"""Tests for the PRIDE import: API parsing, the download job and SDRF settings."""
import functools
import http.server
import os
import stat
import sys
import threading
from pathlib import Path

import pytest

from src import pride

SDRF_HEADER = [
    "source name", "characteristics[organism]", "comment[data file]", "comment[label]",
    "comment[cleavage agent details]", "comment[modification parameters]",
    "comment[modification parameters]", "comment[modification parameters]",
    "comment[precursor mass tolerance]", "comment[fragment mass tolerance]",
    "factor value[disease]",
]
MODS = [
    "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed",
    "NT=Oxidation;MT=Variable;TA=M;AC=UNIMOD:35",
    "NT=Acetyl;AC=UNIMOD:1;PP=Protein N-term;MT=variable",
]


def write_sdrf(path, rows):
    lines = ["\t".join(SDRF_HEADER)] + ["\t".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def lfq_row(data_file, group, frag="20 ppm"):
    return ["s", "Homo sapiens", data_file, "AC=MS:1002038;NT=label free sample",
            "AC=MS:1001251;NT=Trypsin", *MODS, "10 ppm", frag, group]


@pytest.mark.parametrize("text, expected", [
    ("PXD000001", "PXD000001"),
    (" pxd012345 ", "PXD012345"),
    ("https://www.ebi.ac.uk/pride/archive/projects/PXD044981", "PXD044981"),
])
def test_normalize_accession(text, expected):
    assert pride.normalize_accession(text) == expected


def test_normalize_accession_rejects_other_text():
    with pytest.raises(ValueError):
        pride.normalize_accession("MSV000012")


@pytest.mark.parametrize("name, kind", [
    ("run1.raw", pride.KIND_RAW),
    ("run1.RAW", pride.KIND_RAW),
    ("run1.mzML.gz", pride.KIND_MZML),
    ("human.fasta", pride.KIND_FASTA),
    ("PXD000001.sdrf.tsv", pride.KIND_SDRF),
    ("txt.zip", pride.KIND_OTHER),
])
def test_file_kind(name, kind):
    assert pride.file_kind(name) == kind


@pytest.mark.parametrize("value, expected", [
    ("20 ppm", (20.0, "ppm")),
    ("0.02 Da", (0.02, "Da")),
    ("10ppm", (10.0, "ppm")),
    ("20 mmu", (0.02, "Da")),
    ("not available", None),
])
def test_parse_tolerance(value, expected):
    assert pride.parse_tolerance(value) == expected


def test_modification_names_without_validation():
    assert pride.sdrf_modification_names("NT=Phospho;TA=S,T,Y;MT=variable") == [
        "Phospho (S)", "Phospho (T)", "Phospho (Y)"]
    assert pride.sdrf_modification_names("NT=Acetyl;PP=Protein N-term") == ["Acetyl (Protein N-term)"]
    assert pride.sdrf_modification_names("NT=Gln->pyro-Glu;PP=Any N-term;TA=Q") == [
        "Gln->pyro-Glu (N-term Q)"]


def test_sdrf_settings_lfq(tmp_path):
    pytest.importorskip("pyopenms")
    sdrf = pride.read_sdrf(write_sdrf(tmp_path / "a.sdrf.tsv", [
        lfq_row("A_1.raw", "healthy"), lfq_row("A_2.raw", "tumor"), lfq_row("A_3.raw", "tumor"),
    ]))
    settings = pride.sdrf_settings(sdrf, ["A_1.mzML", "A_2.mzML"])
    assert settings["mode"] == "LFQ"
    assert settings["groups"] == {"A_1.mzML": "healthy", "A_2.mzML": "tumor"}
    comet = settings["comet"]
    assert comet["enzyme"] == "Trypsin"
    assert comet["fixed_modifications"] == "Carbamidomethyl (C)"
    assert comet["variable_modifications"] == "Oxidation (M)\nAcetyl (Protein N-term)"
    assert (comet["precursor_mass_tolerance"], comet["precursor_error_units"]) == (10.0, "ppm")
    # Comet bins fragments in Da: 20 ppm becomes a high-resolution 0.02 Da search
    assert (comet["fragment_mass_tolerance"], comet["instrument"], comet["fragment_bin_offset"]) == (
        0.02, "high_res", 0.0)
    assert any("no imported mzML" in n for n in settings["notes"])


def test_sdrf_settings_low_res_fragments(tmp_path):
    sdrf = pride.read_sdrf(write_sdrf(tmp_path / "a.sdrf.tsv", [lfq_row("A_1.raw", "x", "0.5 Da")]))
    comet = pride.sdrf_settings(sdrf)["comet"]
    assert (comet["fragment_mass_tolerance"], comet["instrument"], comet["fragment_bin_offset"]) == (
        0.5, "low_res", 0.4)


def test_sdrf_settings_tmt(tmp_path):
    rows = []
    channels = pride.TMT_CHANNELS["tmt10plex"]
    for i, channel in enumerate(channels):
        row = lfq_row("T_1.raw", "ctrl" if i < 5 else "treated")
        row[3] = f"TMT{channel}"
        rows.append(row)
    settings = pride.sdrf_settings(pride.read_sdrf(write_sdrf(tmp_path / "t.sdrf.tsv", rows)))
    assert settings["mode"] == "TMT"
    assert settings["tmt_type"] == "tmt10plex"
    assert settings["tmt_groups"] == {i: ("ctrl" if i < 5 else "treated") for i in range(10)}


def test_read_sdrf_needs_data_file_column(tmp_path):
    path = tmp_path / "bad.tsv"
    path.write_text("source name\tcharacteristics[organism]\na\tb\n")
    with pytest.raises(ValueError):
        pride.read_sdrf(path)


class FakeResponse:
    def __init__(self, data, status=200):
        self.data, self.status_code = data, status

    def json(self):
        return self.data

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(self.status_code)


class FakeSession:
    def __init__(self, responses):
        self.responses, self.calls = responses, []

    def get(self, url, params=None, timeout=None):
        self.calls.append((url, params))
        return self.responses.pop(0)


def file_entry(name, size=10, category="RAW"):
    return {
        "fileName": name,
        "fileSizeBytes": size,
        "fileCategory": {"value": category},
        "publicFileLocations": [
            {"name": "Aspera Protocol", "value": f"prd_ascp@fasp.ebi.ac.uk:pride/{name}"},
            {"name": "FTP Protocol", "value": f"ftp://ftp.pride.ebi.ac.uk/pride/data/archive/{name}"},
        ],
    }


def test_fetch_files_pages_and_uses_https(monkeypatch):
    monkeypatch.setattr(pride, "FILES_PAGE_SIZE", 2)
    session = FakeSession([
        FakeResponse([file_entry("a.raw"), file_entry("b.raw")]),
        FakeResponse([file_entry("x.sdrf.tsv", category="EXPERIMENTAL DESIGN")]),
    ])
    files = pride.fetch_files("PXD000001", session)
    assert [f.name for f in files] == ["a.raw", "b.raw", "x.sdrf.tsv"]
    assert files[0].url == "https://ftp.pride.ebi.ac.uk/pride/data/archive/a.raw"
    assert files[2].kind == pride.KIND_SDRF
    assert [c[1]["page"] for c in session.calls] == [0, 1]


def test_fetch_project(monkeypatch):
    session = FakeSession([FakeResponse({
        "accession": "PXD000001", "title": "T",
        "organisms": [{"name": "Homo sapiens (human)"}],
        "instruments": [{"name": "Q Exactive"}, {"name": "Q Exactive"}],
        "references": [{"referenceLine": "Paper", "doi": "10.1/x"}],
    })])
    project = pride.fetch_project("PXD000001", session)
    assert project["organisms"] == ["Homo sapiens (human)"]
    assert project["instruments"] == ["Q Exactive"]
    assert project["references"] == ["Paper doi:10.1/x"]


def test_fetch_project_not_found():
    with pytest.raises(LookupError):
        pride.fetch_project("PXD999999", FakeSession([FakeResponse({}, 404)]))


@pytest.fixture
def file_server(tmp_path):
    root = tmp_path / "served"
    root.mkdir()
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(root))
    handler.log_message = lambda *a: None
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    yield root, f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()


@pytest.fixture
def fake_file_converter(tmp_path, monkeypatch):
    """A FileConverter on PATH that writes a stub mzML."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    script = bin_dir / "FileConverter"
    script.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        "args = sys.argv[1:]\n"
        "out = args[args.index('-out') + 1]\n"
        "open(out, 'w').write('<mzML/>')\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")


def test_run_import(tmp_path, file_server, fake_file_converter):
    root, base = file_server
    (root / "run1.mzML").write_text("<mzML>1</mzML>")
    (root / "run2.raw").write_bytes(b"raw bytes")
    (root / "db.fasta").write_text(">P1\nPEPTIDE\n")
    write_sdrf(root / "PXD000001.sdrf.tsv", [lfq_row("run1.raw", "a")])

    def entry(name):
        size = (root / name).stat().st_size
        return {"name": name, "size": size, "url": f"{base}/{name}", "kind": pride.file_kind(name)}

    job, mzml_dir, fasta_dir = tmp_path / "job", tmp_path / "mzML-files", tmp_path / "fasta-file"
    job.mkdir()
    mzml_dir.mkdir()
    (mzml_dir / "example.mzML").write_text("example")  # what an empty upload tab shows
    files = [entry(n) for n in ("run1.mzML", "run2.raw", "db.fasta", "PXD000001.sdrf.tsv")]
    files.append({"name": "missing.raw", "size": 3, "url": f"{base}/missing.raw", "kind": pride.KIND_RAW})
    pride.run_import(str(job), "PXD000001", files, str(mzml_dir), str(fasta_dir),
                     example_files={str(mzml_dir): ["example.mzML"]})

    status = pride.read_status(job)
    assert status["state"] == "done"
    assert status["done_files"] == 5
    assert sorted(p.name for p in mzml_dir.iterdir()) == ["run1.mzML", "run2.mzML"]
    assert (fasta_dir / "db.fasta").exists()
    assert status["sdrf"] == [str(job / "sdrf" / "PXD000001.sdrf.tsv")]
    assert len(status["errors"]) == 1 and status["errors"][0].startswith("missing.raw")
    # RAW files are removed after conversion unless asked to keep them
    assert not (job / "downloads" / "PXD000001" / "run2.raw").exists()


def test_download_resumes_partial_file(tmp_path, file_server):
    root, base = file_server
    (root / "big.mzML").write_bytes(b"0123456789")
    dest = tmp_path / "big.mzML"
    (tmp_path / "big.mzML.part").write_bytes(b"01234")
    # SimpleHTTPRequestHandler ignores Range and answers 200: the download restarts
    pride.download(f"{base}/big.mzML", dest, 10)
    assert dest.read_bytes() == b"0123456789"
    assert not (tmp_path / "big.mzML.part").exists()


def test_download_rejects_wrong_size(tmp_path, file_server):
    root, base = file_server
    (root / "a.mzML").write_bytes(b"short")
    with pytest.raises(IOError):
        pride.download(f"{base}/a.mzML", tmp_path / "a.mzML", 999)
