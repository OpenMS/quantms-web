"""Import PRIDE Archive datasets into the workflow's input files.

Everything here runs without Streamlit: the download job is started as its own
process by ``src/pride_ui.py`` and reports progress through a status file, so
the page can be left and reopened while a dataset is downloading.

PRIDE Archive REST API (v3): https://www.ebi.ac.uk/pride/ws/archive/v3/webjars/swagger-ui/index.html
SDRF-Proteomics: https://github.com/bigbio/proteomics-sample-metadata
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import requests

PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
PRIDE_PROJECT_PAGE = "https://www.ebi.ac.uk/pride/archive/projects/{accession}"
REQUEST_TIMEOUT = 60
FILES_PAGE_SIZE = 100
CHUNK_SIZE = 1 << 20

ACCESSION_RE = re.compile(r"^P[XR]D\d{6}$")

# File kinds the import knows what to do with.
KIND_MZML = "mzML"
KIND_RAW = "Thermo RAW"
KIND_FASTA = "FASTA"
KIND_SDRF = "SDRF"
KIND_OTHER = "other"


def normalize_accession(text: str) -> str:
    """Return a PRIDE accession (PXD/PRD + 6 digits) from user input.

    Accepts the bare accession, lower case, or a PRIDE project URL.
    """
    match = re.search(r"P[XR]D\d{6}", (text or "").strip().upper())
    if not match:
        raise ValueError(f"'{text}' is not a PRIDE accession (e.g. PXD000001).")
    return match.group(0)


@dataclass
class PrideFile:
    name: str
    category: str
    size: int
    url: str
    kind: str

    @property
    def size_mb(self) -> float:
        return round(self.size / 1e6, 1)


def file_kind(name: str, category: str = "") -> str:
    """Classify a dataset file by name (and PRIDE category for SDRF files)."""
    lower = name.lower()
    if lower.endswith(".gz"):
        lower = lower[:-3]
    if lower.endswith(".mzml"):
        return KIND_MZML
    if lower.endswith(".raw"):
        return KIND_RAW
    if lower.endswith((".fasta", ".fa", ".faa")):
        return KIND_FASTA
    if "sdrf" in lower and lower.endswith((".tsv", ".txt", ".csv")):
        return KIND_SDRF
    if "EXPERIMENTAL DESIGN" in category.upper() and lower.endswith(".tsv"):
        return KIND_SDRF
    return KIND_OTHER


def _session() -> requests.Session:
    session = requests.Session()
    session.headers["Accept"] = "application/json"
    return session


def _names(entries) -> list[str]:
    """Names from PRIDE's lists of CV params / plain strings."""
    names = []
    for entry in entries or []:
        name = entry.get("name") if isinstance(entry, dict) else entry
        if name and name not in names:
            names.append(str(name))
    return names


def fetch_project(accession: str, session: requests.Session = None) -> dict:
    """Project metadata, reduced to what the import page shows."""
    session = session or _session()
    response = session.get(f"{PRIDE_API}/projects/{accession}", timeout=REQUEST_TIMEOUT)
    if response.status_code == 404:
        raise LookupError(f"{accession} was not found in PRIDE (or is not public yet).")
    response.raise_for_status()
    data = response.json()
    references = []
    for ref in data.get("references") or []:
        line = ref.get("referenceLine") or ""
        doi = ref.get("doi") or ""
        if doi and doi not in line:
            line = f"{line} doi:{doi}".strip()
        if line:
            references.append(line)
    return {
        "accession": data.get("accession", accession),
        "title": data.get("title", ""),
        "description": data.get("projectDescription", ""),
        "sample_protocol": data.get("sampleProcessingProtocol", ""),
        "data_protocol": data.get("dataProcessingProtocol", ""),
        "submission_date": data.get("submissionDate", ""),
        "publication_date": data.get("publicationDate", ""),
        "organisms": _names(data.get("organisms")),
        "instruments": _names(data.get("instruments")),
        "ptms": _names(data.get("identifiedPTMStrings")),
        "quantification": _names(data.get("quantificationMethods")),
        "experiment_types": _names(data.get("experimentTypes")),
        "keywords": _names(data.get("keywords")),
        "references": references,
        "url": PRIDE_PROJECT_PAGE.format(accession=accession),
    }


def _download_url(entry: dict) -> str:
    """HTTPS download link for a PRIDE file entry.

    PRIDE lists FTP and Aspera locations; the FTP server serves the same paths
    over HTTPS, which works through web proxies and needs no FTP client.
    """
    locations = entry.get("publicFileLocations") or []
    for location in locations:
        value = location.get("value", "")
        if value.startswith("ftp://"):
            return "https://" + value[len("ftp://"):]
        if value.startswith("https://"):
            return value
    return ""


def fetch_files(accession: str, session: requests.Session = None) -> list[PrideFile]:
    """All files of a project, with HTTPS download links."""
    session = session or _session()
    files, page = [], 0
    while True:
        response = session.get(
            f"{PRIDE_API}/projects/{accession}/files",
            params={"pageSize": FILES_PAGE_SIZE, "page": page},
            timeout=REQUEST_TIMEOUT,
        )
        response.raise_for_status()
        data = response.json()
        # v3 returns a plain list; older HAL responses nest it in _embedded
        entries = data if isinstance(data, list) else data.get("_embedded", {}).get("files", [])
        for entry in entries:
            name = entry.get("fileName", "")
            category = entry.get("fileCategory") or {}
            category = category.get("value", "") if isinstance(category, dict) else str(category)
            files.append(PrideFile(
                name=name,
                category=category,
                size=int(entry.get("fileSizeBytes") or 0),
                url=_download_url(entry),
                kind=file_kind(name, category),
            ))
        if len(entries) < FILES_PAGE_SIZE:
            break
        page += 1
    return files


def raw_conversion_available() -> tuple[bool, str]:
    """Whether FileConverter can turn Thermo RAW files into mzML here.

    FileConverter calls the ThermoRawFileParser shipped with OpenMS, which on
    Linux and macOS runs on Mono.
    """
    if not shutil.which("FileConverter"):
        return False, "The OpenMS FileConverter tool is not installed."
    if sys.platform != "win32" and not shutil.which("mono"):
        return False, "Mono is not installed, which ThermoRawFileParser needs on Linux and macOS."
    return True, ""


# ---------------------------------------------------------------------------
# Download job
# ---------------------------------------------------------------------------

def status_path(job_dir: Path) -> Path:
    return Path(job_dir, "status.json")


def read_status(job_dir: Path) -> dict:
    try:
        with open(status_path(job_dir), "r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def write_status(job_dir: Path, status: dict) -> None:
    """Write the status file atomically, so the page never reads half a file."""
    path = status_path(job_dir)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(status, f, indent=2)
    os.replace(tmp, path)


def is_running(job_dir: Path) -> bool:
    status = read_status(job_dir)
    if status.get("state") != "running":
        return False
    try:
        os.kill(int(status.get("pid", 0)), 0)
    except (OSError, ValueError):
        return False
    return True


def cancel(job_dir: Path) -> None:
    status = read_status(job_dir)
    try:
        os.kill(int(status.get("pid", 0)), signal.SIGTERM)
    except (OSError, ValueError):
        pass


def _log(job_dir: Path, message: str) -> None:
    with open(Path(job_dir, "import.log"), "a", encoding="utf-8") as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {message}\n")


def download(url: str, dest: Path, expected_size: int = 0, on_progress=None,
             session: requests.Session = None) -> Path:
    """Stream ``url`` to ``dest``, resuming a partial ``.part`` file if there is one."""
    if dest.exists() and (not expected_size or dest.stat().st_size == expected_size):
        if on_progress:
            on_progress(dest.stat().st_size)
        return dest
    session = session or requests.Session()
    part = dest.with_name(dest.name + ".part")
    offset = part.stat().st_size if part.exists() else 0
    headers = {"Range": f"bytes={offset}-"} if offset else {}
    with session.get(url, stream=True, headers=headers, timeout=REQUEST_TIMEOUT) as response:
        if response.status_code == 416:  # the .part file is already complete
            response.close()
        else:
            response.raise_for_status()
            if offset and response.status_code != 206:
                offset = 0  # server ignored the range: start over
            with open(part, "ab" if offset else "wb") as f:
                done = offset
                for chunk in response.iter_content(CHUNK_SIZE):
                    f.write(chunk)
                    done += len(chunk)
                    if on_progress:
                        on_progress(done)
    if expected_size and part.stat().st_size != expected_size:
        raise IOError(
            f"{dest.name}: downloaded {part.stat().st_size} bytes, PRIDE lists {expected_size}."
        )
    os.replace(part, dest)
    return dest


def _gunzip(path: Path) -> Path:
    out = path.with_suffix("")
    with gzip.open(path, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst, CHUNK_SIZE)
    path.unlink()
    return out


def run_import(job_dir: str, accession: str, files: list[dict], mzml_dir: str,
               fasta_dir: str, keep_raw: bool = False, example_files: dict = None) -> None:
    """Download ``files`` and place them in the workflow's input folders.

    mzML files go to ``mzml_dir``; Thermo RAW files are converted there with
    FileConverter; FASTA files go to ``fasta_dir``; SDRF files stay in
    ``job_dir`` for the page to apply. Runs as its own process.

    ``example_files`` maps an input folder to the example file names the upload
    tab puts into it while it is empty; if a folder holds just those, they are
    removed before the first imported file goes in.
    """
    example_files = example_files or {}
    checked = set()

    def target_in(directory: Path, name: str) -> Path:
        """``directory/name``; replaces example data and earlier imports of the file."""
        if directory not in checked:
            checked.add(directory)
            examples = sorted(example_files.get(str(directory), []))
            if examples and sorted(p.name for p in directory.iterdir()) == examples:
                for example in examples:
                    Path(directory, example).unlink()
        target = Path(directory, name)
        if target.exists():
            target.unlink()
        return target

    job_dir, mzml_dir, fasta_dir = Path(job_dir), Path(mzml_dir), Path(fasta_dir)
    downloads = Path(job_dir, "downloads", accession)
    downloads.mkdir(parents=True, exist_ok=True)
    mzml_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)

    status = {
        "state": "running",
        "pid": os.getpid(),
        "accession": accession,
        "started": time.time(),
        "total_files": len(files),
        "done_files": 0,
        "bytes_total": sum(int(f.get("size") or 0) for f in files),
        "bytes_done": 0,
        "current": "",
        "step": "",
        "imported": [],
        "sdrf": [],
        "errors": [],
    }
    child = {"process": None}

    def on_term(signum, frame):
        if child["process"] is not None:
            child["process"].terminate()
        status["state"] = "cancelled"
        status["finished"] = time.time()
        write_status(job_dir, status)
        _log(job_dir, "Cancelled.")
        os._exit(0)

    signal.signal(signal.SIGTERM, on_term)
    write_status(job_dir, status)
    _log(job_dir, f"Importing {len(files)} file(s) from {accession}.")

    session = requests.Session()
    bytes_before = 0
    last_write = [0.0]

    def on_progress(done):
        status["bytes_done"] = bytes_before + done
        if time.time() - last_write[0] > 1:
            write_status(job_dir, status)
            last_write[0] = time.time()

    for entry in files:
        name, kind = entry["name"], entry["kind"]
        status["current"], status["step"] = name, "download"
        write_status(job_dir, status)
        try:
            path = download(entry["url"], Path(downloads, name), int(entry.get("size") or 0),
                            on_progress, session)
            bytes_before += int(entry.get("size") or path.stat().st_size)
            if path.suffix.lower() == ".gz":
                status["step"] = "decompress"
                write_status(job_dir, status)
                path = _gunzip(path)

            if kind == KIND_MZML:
                target = target_in(mzml_dir, path.name)
                shutil.move(str(path), target)
            elif kind == KIND_RAW:
                status["step"] = "convert"
                write_status(job_dir, status)
                target = target_in(mzml_dir, path.stem + ".mzML")
                command = ["FileConverter", "-in", str(path), "-out", str(target)]
                _log(job_dir, "Running " + " ".join(command))
                child["process"] = subprocess.Popen(
                    command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
                )
                output, _ = child["process"].communicate()
                returncode = child["process"].returncode
                child["process"] = None
                _log(job_dir, output.strip())
                if returncode != 0 or not target.exists():
                    target.unlink(missing_ok=True)
                    raise RuntimeError(f"FileConverter failed (exit code {returncode}).")
                if not keep_raw:
                    path.unlink()
            elif kind == KIND_FASTA:
                target = target_in(fasta_dir, path.name)
                shutil.move(str(path), target)
            elif kind == KIND_SDRF:
                target = Path(job_dir, "sdrf", path.name)
                target.parent.mkdir(exist_ok=True)
                shutil.move(str(path), target)
                status["sdrf"].append(str(target))
            else:
                target = path
            status["imported"].append(str(target))
            _log(job_dir, f"{name}: imported as {target}")
        except Exception as e:  # keep going with the remaining files
            status["errors"].append(f"{name}: {e}")
            _log(job_dir, f"{name}: ERROR {e}")
        status["done_files"] += 1
        status["bytes_done"] = bytes_before
        write_status(job_dir, status)

    status["state"] = "failed" if status["errors"] and not status["imported"] else "done"
    status["current"], status["step"] = "", ""
    status["finished"] = time.time()
    write_status(job_dir, status)
    _log(job_dir, f"Finished: {len(status['imported'])} imported, {len(status['errors'])} failed.")


# ---------------------------------------------------------------------------
# SDRF
# ---------------------------------------------------------------------------

@dataclass
class Sdrf:
    columns: list[str]
    rows: list[list[str]]

    def values(self, column: str) -> list[list[str]]:
        """Per row, every non-empty value of the (possibly repeated) column."""
        idx = [i for i, c in enumerate(self.columns) if c == column]
        return [[row[i].strip() for i in idx if i < len(row) and row[i].strip()] for row in self.rows]

    def column(self, column: str) -> list[str]:
        """First value of ``column`` per row ('' if absent)."""
        return [(v[0] if v else "") for v in self.values(column)]

    def has(self, column: str) -> bool:
        return column in self.columns


def read_sdrf(path) -> Sdrf:
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        lines = [row for row in csv.reader(f, delimiter="\t") if any(c.strip() for c in row)]
    if not lines:
        raise ValueError("The SDRF file is empty.")
    columns = [re.sub(r"\s+", " ", c.strip().lower()) for c in lines[0]]
    if "comment[data file]" not in columns:
        raise ValueError("The SDRF file has no 'comment[data file]' column.")
    return Sdrf(columns, lines[1:])


def parse_key_values(value: str) -> dict:
    """'NT=Oxidation;AC=UNIMOD:35;TA=M' -> {'NT': 'Oxidation', ...}. A bare value becomes NT."""
    result = {}
    for part in value.split(";"):
        if "=" in part:
            key, val = part.split("=", 1)
            result[key.strip().upper()] = val.strip()
        elif part.strip():
            result.setdefault("NT", part.strip())
    return result


def parse_tolerance(value: str) -> tuple[float, str] | None:
    """'20 ppm' -> (20.0, 'ppm'); '0.02 Da' -> (0.02, 'Da')."""
    match = re.match(r"^\s*([0-9.]+)\s*(ppm|da|th|mmu)?\s*$", value or "", re.IGNORECASE)
    if not match:
        return None
    unit = (match.group(2) or "Da").lower()
    amount = float(match.group(1))
    if unit == "mmu":
        return amount / 1000.0, "Da"
    return amount, "ppm" if unit == "ppm" else "Da"


def _openms_mods():
    """Modification and enzyme names known to OpenMS, or None without pyOpenMS."""
    try:
        import pyopenms as oms
    except ImportError:
        return None, None, {}
    mods = []
    oms.ModificationsDB().getAllSearchModifications(mods)
    mods = [m.decode() if isinstance(m, bytes) else m for m in mods]
    enzymes = []
    oms.ProteaseDB().getAllNames(enzymes)
    enzymes = [e.decode() if isinstance(e, bytes) else e for e in enzymes]
    # UniMod accession -> OpenMS modification id, to resolve SDRF names that differ
    by_unimod = {}
    db = oms.ModificationsDB()
    for full_id in mods:
        try:
            mod = db.getModification(full_id)
            unimod = mod.getUniModAccession()
            unimod = unimod.decode() if isinstance(unimod, bytes) else unimod
            mod_id = mod.getId()
            by_unimod.setdefault(unimod.upper(), mod_id.decode() if isinstance(mod_id, bytes) else mod_id)
        except Exception:
            continue
    return set(mods), enzymes, by_unimod


_TERMINI = {
    "protein n-term": "Protein N-term",
    "any n-term": "N-term",
    "protein c-term": "Protein C-term",
    "any c-term": "C-term",
}


def sdrf_modification_names(value: str, known=None, by_unimod=None) -> list[str]:
    """OpenMS/Comet names for one SDRF modification entry, e.g. 'Oxidation (M)'."""
    kv = parse_key_values(value)
    name = kv.get("NT", "")
    accession = kv.get("AC", "").upper().replace(" ", "")
    if known is not None and by_unimod and accession in by_unimod:
        name = by_unimod[accession]
    if not name:
        return []
    position = _TERMINI.get(kv.get("PP", "").strip().lower(), "")
    sites = [s.strip() for s in kv.get("TA", "").split(",") if s.strip() and s.strip().upper() != "X"]
    if position:
        candidates = [f"{name} ({position} {s})" for s in sites] or [f"{name} ({position})"]
    else:
        candidates = [f"{name} ({s})" for s in sites]
    if known is None:
        return candidates
    resolved = []
    for candidate in candidates:
        match = next((m for m in known if m.lower() == candidate.lower()), None)
        if match:
            resolved.append(match)
    return resolved


def _modification_type(value: str) -> str:
    return parse_key_values(value).get("MT", "variable").strip().lower()


def _most_common(values: list[str]) -> str:
    values = [v for v in values if v]
    return max(set(values), key=values.count) if values else ""


# Reporter channel order per plex, as IsobaricAnalyzer numbers them.
TMT_CHANNELS = {
    "tmt6plex": ["126", "127", "128", "129", "130", "131"],
    "tmt10plex": ["126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131"],
    "tmt11plex": ["126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C"],
    "tmt16plex": ["126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N",
                  "131C", "132N", "132C", "133N", "133C", "134N"],
    "tmt18plex": ["126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N",
                  "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N"],
}


def _tmt_channel(label: str) -> str:
    """'TMT127N' -> '127N'; 'TMT131' -> '131'."""
    match = re.search(r"(\d{3}[NC]?)", label.upper())
    return match.group(1) if match else ""


def sdrf_settings(sdrf: Sdrf, mzml_names: list[str] = ()) -> dict:
    """Workflow settings an SDRF file implies.

    Returns a dict with ``mode`` ('LFQ'/'TMT'), ``groups`` (mzML file name ->
    group, LFQ), ``tmt_type`` and ``tmt_groups`` (channel index -> group, TMT),
    ``comet`` (CometAdapter parameter -> value) and ``notes`` (what could not be
    applied).
    """
    known, enzymes, by_unimod = _openms_mods()
    notes = []

    # Group = factor values of the run (several factors joined with "_")
    factor_columns = [c for c in dict.fromkeys(sdrf.columns) if c.startswith("factor value[")]
    if not factor_columns:
        notes.append("The SDRF has no factor value columns, so no sample groups are set.")
    groups_per_row = [
        "_".join(v for c in factor_columns for v in sdrf.values(c)[i][:1])
        for i in range(len(sdrf.rows))
    ]
    data_files = sdrf.column("comment[data file]")
    labels = sdrf.column("comment[label]")
    tmt_rows = [i for i, label in enumerate(labels) if "TMT" in label.upper()]
    mode = "TMT" if tmt_rows else "LFQ"

    settings = {"mode": mode, "groups": {}, "tmt_type": "", "tmt_groups": {}, "comet": {},
                "notes": notes}

    if mode == "LFQ":
        by_stem = {Path(n).stem.lower(): n for n in mzml_names}
        unmatched = []
        for data_file, group in zip(data_files, groups_per_row):
            stem = Path(data_file).name.lower()
            for ext in (".gz", ".raw", ".mzml"):
                stem = stem[: -len(ext)] if stem.endswith(ext) else stem
            target = by_stem.get(stem, f"{Path(data_file).name.split('.')[0]}.mzML")
            if mzml_names and stem not in by_stem:
                unmatched.append(data_file)
                continue
            if group:
                settings["groups"][target] = group
        if unmatched:
            notes.append(f"{len(unmatched)} SDRF data file(s) have no imported mzML file and were skipped.")
    else:
        plex = len({_tmt_channel(labels[i]) for i in tmt_rows})
        tmt_type = next((t for t, ch in TMT_CHANNELS.items() if len(ch) >= plex), "tmt18plex")
        if plex in (6, 10, 11, 16, 18):
            tmt_type = f"tmt{plex}plex"
        settings["tmt_type"] = tmt_type
        channels = TMT_CHANNELS[tmt_type]
        per_channel = {}
        for i in tmt_rows:
            channel = _tmt_channel(labels[i])
            if channel in channels:
                per_channel.setdefault(channels.index(channel), set()).add(groups_per_row[i])
        for index, values in sorted(per_channel.items()):
            if len(values) == 1 and next(iter(values)):
                settings["tmt_groups"][index] = next(iter(values))
        if any(len(v) > 1 for v in per_channel.values()):
            notes.append("Channels carry different groups in different runs; only consistent channels were set.")

    # Search settings: the value most runs share
    comet = settings["comet"]
    enzyme_values = [parse_key_values(v[0]).get("NT", "") if v else "" for v in sdrf.values("comment[cleavage agent details]")]
    enzyme = _most_common(enzyme_values)
    if enzyme:
        if enzymes is None:
            comet["enzyme"] = enzyme
        else:
            match = next((e for e in enzymes if e.lower() == enzyme.lower()), None)
            if match:
                comet["enzyme"] = match
            else:
                notes.append(f"Enzyme '{enzyme}' is not known to OpenMS and was not set.")

    fixed, variable = [], []
    for entries in sdrf.values("comment[modification parameters]"):
        for entry in entries:
            names = sdrf_modification_names(entry, known, by_unimod)
            if not names:
                notes.append(f"Modification '{entry}' could not be matched to an OpenMS modification.")
            target = fixed if _modification_type(entry) == "fixed" else variable
            for name in names:
                if name not in target:
                    target.append(name)
    if fixed or variable:
        comet["fixed_modifications"] = "\n".join(fixed)
        comet["variable_modifications"] = "\n".join(variable)

    precursor = parse_tolerance(_most_common(sdrf.column("comment[precursor mass tolerance]")))
    if precursor:
        comet["precursor_mass_tolerance"], comet["precursor_error_units"] = precursor
    fragment = parse_tolerance(_most_common(sdrf.column("comment[fragment mass tolerance]")))
    if fragment:
        amount, unit = fragment
        # Comet bins fragments in Da; ppm tolerances mean a high-resolution MS2
        if unit == "ppm":
            amount = max(round(amount * 1e-3, 3), 0.01)
        high_res = unit == "ppm" or amount <= 0.1
        comet["fragment_mass_tolerance"] = amount
        comet["instrument"] = "high_res" if high_res else "low_res"
        comet["fragment_bin_offset"] = 0.0 if high_res else 0.4

    # Notes repeat once per run; keep each once
    settings["notes"] = list(dict.fromkeys(notes))
    return settings


def settings_preview(settings: dict) -> list[dict]:
    """Rows describing what applying ``settings`` changes, for display."""
    rows = [{"Setting": "Analysis mode", "Value": settings["mode"]}]
    for name, group in settings["groups"].items():
        rows.append({"Setting": f"Group for {name}", "Value": group})
    if settings["tmt_type"]:
        rows.append({"Setting": "Isobaric type", "Value": settings["tmt_type"]})
    for index, group in settings["tmt_groups"].items():
        rows.append({"Setting": f"Group for channel {index + 1}", "Value": group})
    for key, value in settings["comet"].items():
        value = ", ".join(value.split("\n")) if isinstance(value, str) else value
        rows.append({"Setting": f"Comet {key}", "Value": str(value) or "(none)"})
    return rows


def as_dicts(files: list[PrideFile]) -> list[dict]:
    return [asdict(f) for f in files]
