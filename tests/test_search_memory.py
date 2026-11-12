"""The memory estimate of the database search (src/search_memory.py)."""
from types import SimpleNamespace

import pytest

from src import search_memory as sm

GB = sm.GB
MB = sm.MB


def test_fasta_counts(tmp_path):
    fasta = tmp_path / "db.fasta"
    fasta.write_text(">a\nPEPTIDE\nKR\n>b\nACDE\n")
    assert sm.fasta_counts(fasta) == (2, 13)


@pytest.mark.parametrize("variable, n", [
    ("Oxidation (M)\nAcetyl (Protein N-term)", 1),
    (["Oxidation (M)", "Phospho (S)", "Phospho (T)"], 3),
    ("", 0), (None, 0),
])
def test_index_modifications(variable, n):
    assert sm.index_modifications(variable) == n


def test_specificity():
    assert sm.specificity({}) == "full"
    assert sm.specificity({"Search:peptide:enzyme_specificity": "semi"}) == "semi"
    assert sm.specificity({"Search:enzyme": "unspecific cleavage"}) == "none"


def test_index_grows_with_proteins_specificity_and_mods():
    human = 20_000
    tryptic = sm.index_bytes(human, human * 555, "full", 1)
    assert tryptic == pytest.approx(0.5 * GB)
    assert sm.index_bytes(human, 0, "semi", 1) == pytest.approx(15 * tryptic)
    assert sm.index_bytes(human, 0, "full", 3) == pytest.approx(4 * tryptic)
    # A chunk costs twice the unchunked index per protein, at most the whole
    assert sm.index_bytes(human, 0, "full", 1, chunk=1000) < tryptic
    assert sm.index_bytes(human, 0, "full", 1, chunk=19_999) == pytest.approx(tryptic)
    nonspecific = sm.index_bytes(4000, 1_460_000, "none", 1)
    assert nonspecific == pytest.approx(9.6 * 1.46 * GB)
    assert sm.index_bytes(4000, 1_460_000, "none", 1, chunk=250) < nonspecific


def test_growth_levels_off():
    assert sm.growth(1, 50 * MB) == 0
    assert sm.growth(5, 50 * MB) == pytest.approx(4 * sm.GROWTH_SMALL)
    assert sm.growth(9, 1000 * MB) == pytest.approx(4 * sm.GROWTH_LARGE + sm.GROWTH_LARGE)


def test_processes_side_by_side_each_hold_an_index():
    values = {"Search:peptide:enzyme_specificity": "semi", "Search:modifications:variable": "Oxidation (M)"}
    sizes = [500 * MB] * 8
    one = sm.search_memory(sizes, [list(range(8))], 30 * MB, 20_000, 0, values, 16)
    four = sm.search_memory(sizes, [[0, 1], [2, 3], [4, 5], [6, 7]], 30 * MB, 20_000, 0, values, 64)
    assert one["cause"] == "database"
    assert four["bytes"] > 3.5 * one["bytes"]
    # A semi-specific human search is about 7.5 GB per process
    assert 7.5 * GB < one["bytes"] < 10 * GB
    small = sm.search_memory(sizes, [list(range(8))], 30 * MB, 4000, 0, {}, 16)
    assert small["cause"] == "files"


def _workflow(tmp_path, threads=64):
    from src.WorkflowTest import WorkflowTest

    wf = object.__new__(WorkflowTest)
    wf.executor = SimpleNamespace(_get_max_threads=lambda: threads)
    fasta = tmp_path / "human.fasta"
    fasta.write_text("".join(f">sp|P{i}|X\n{'A' * 555}\n" for i in range(20_000)))
    files = []
    for i in range(8):
        f = tmp_path / f"s{i}.mzML"
        with open(f, "wb") as h:
            h.truncate(400 * MB)  # sparse
        files.append(str(f))
    return wf, {"mzML-files": files, "fasta-file": str(fasta)}


def test_workflow_estimate_and_advice(tmp_path, monkeypatch):
    wf, params = _workflow(tmp_path)
    tryptic = wf.memory_estimate(params)
    assert tryptic and tryptic > 4 * GB  # 4 processes of 16 threads
    params["ProSE"] = {"Search:peptide:enzyme_specificity": "semi"}
    semi = wf.memory_estimate(params)
    assert semi > 4 * tryptic
    assert "chunk size" in wf.memory_advice(params)
    monkeypatch.setattr("src.workflow.settings_io.load_settings", lambda *a: {"online_deployment": True})
    wf.memory_estimate({**params, "ProSE": {}, "fasta-file": params["fasta-file"]})
    wf._memory_cause = "files"
    assert "two or more runs" in wf.memory_advice(params)


def test_tmt_reads_its_own_search_settings(tmp_path):
    wf, params = _workflow(tmp_path)
    params["analysis-mode"] = "TMT"
    params["ProSE"] = {"Search:peptide:enzyme_specificity": "semi"}
    base = wf.memory_estimate(params)
    params["ProSE-TMT"] = {"Search:peptide:enzyme_specificity": "semi"}
    assert wf.memory_estimate(params) > 4 * base


def test_no_estimate_without_files_or_database(tmp_path):
    wf, params = _workflow(tmp_path)
    assert wf.memory_estimate({**params, "mzML-files": []}) is None
    assert wf.memory_estimate({**params, "fasta-file": str(tmp_path / "missing.fasta")}) is None


def test_raw_counts_twice_and_d_folders_are_summed(tmp_path):
    from src.WorkflowTest import _input_bytes

    raw = tmp_path / "a.raw"
    raw.write_bytes(b"x" * 100)
    d = tmp_path / "b.d"
    d.mkdir()
    (d / "analysis.tdf_bin").write_bytes(b"x" * 300)
    assert _input_bytes(raw) == 200
    assert _input_bytes(d) == 300
    assert _input_bytes(tmp_path / "gone.mzML") == 0
