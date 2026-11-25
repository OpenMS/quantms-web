"""A file Percolator cannot rescore is left out; the rest of the study goes on.

execution() runs PercolatorAdapter with allow_failed=True. When it fails on a
file (a blank, a run with too few matches to learn from), that file is dropped
from every list the next steps take, so IDFilter and ProteomicsLFQ get only the
files that were rescored. The tools are faked; the run stops at ProteomicsLFQ,
whose inputs are what the test checks.
"""
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

# Imported inside the test, as in test_search_batches.py.


class FakeExecutor:
    def __init__(self, percolator_fails):
        self.percolator_fails = percolator_fails
        self.failed = []
        self.calls = []

    def convert_to_mzml(self, files):
        return list(files)

    def run_topp(self, tool, io, params=None, tool_instance_name=None, allow_failed=False):
        self.calls.append((tool, io, allow_failed))
        if tool == "PercolatorAdapter":
            self.failed = list(self.percolator_fails)
            if not allow_failed:
                return not self.failed
            return len(self.failed) < len(io["in"])
        self.failed = []
        return tool != "ProteomicsLFQ"  # stop the run there

    def keep_succeeded(self, values):
        return [v for i, v in enumerate(values) if i not in set(self.failed)]


def _workflow(tmp_path, monkeypatch, percolator_fails):
    import src.WorkflowTest as module
    from src.WorkflowTest import WorkflowTest

    names = ["A_1.mzML", "Blank.mzML", "B_1.mzML"]
    mzml = [str(tmp_path / n) for n in names]
    wf = object.__new__(WorkflowTest)
    wf.workflow_dir = tmp_path
    wf.params = {"mzML-files": mzml, "fasta-file": str(tmp_path / "db.fasta"), "analysis-mode": "LFQ"}
    wf.logger = SimpleNamespace(log=lambda *a, **k: None)
    wf.executor = FakeExecutor(percolator_fails)
    wf.file_manager = SimpleNamespace(get_files=lambda files: list(files))
    monkeypatch.setattr(WorkflowTest, "results_dir", property(lambda self: tmp_path / "results"), raising=False)
    monkeypatch.setattr(wf, "prepare_search_database", lambda fasta: (tmp_path / "db.fasta", "DECOY_"), raising=False)
    monkeypatch.setattr(wf, "run_search", lambda *a, **k: True, raising=False)
    monkeypatch.setattr(wf, "search_tolerance", lambda inst: (20.0, True), raising=False)
    monkeypatch.setattr(wf, "_has_identifications", lambda files, stage: True, raising=False)
    monkeypatch.setattr(wf, "add_prediction_features", lambda results, mz, inst: list(results), raising=False)
    monkeypatch.setattr(
        wf, "filter_psms",
        lambda i, o, **k: wf.executor.calls.append(("IDFilter", {"in": i, "out": o}, False)) or True,
        raising=False,
    )
    monkeypatch.setattr(module, "write_spectra_cache", lambda *a, **k: None)
    monkeypatch.setattr(module, "build_psm_views", lambda *a, **k: True)
    monkeypatch.setattr(
        module, "load_design",
        lambda wd, files: pd.DataFrame({"file": [Path(f).name for f in files], "fraction": 1}),
    )
    monkeypatch.setattr(module, "is_fractionated", lambda design: False)
    return wf, names


def _inputs(wf, tool):
    return next(io for t, io, _ in wf.executor.calls if t == tool)


def test_the_file_percolator_fails_on_is_left_out(tmp_path, monkeypatch):
    wf, names = _workflow(tmp_path, monkeypatch, percolator_fails=[1])
    assert wf.execution() is False  # the fake ProteomicsLFQ stops the run

    percolator = next(c for c in wf.executor.calls if c[0] == "PercolatorAdapter")
    assert percolator[2] is True  # allow_failed
    assert len(percolator[1]["in"]) == 3

    idfilter = _inputs(wf, "IDFilter")
    assert [Path(f).name for f in idfilter["in"]] == ["A_1_per.idXML", "B_1_per.idXML"]
    lfq = _inputs(wf, "ProteomicsLFQ")
    assert [Path(f).name for f in lfq["in"][0]] == ["A_1.mzML", "B_1.mzML"]
    assert [Path(f).name for f in lfq["ids"][0]] == ["A_1_filter.idXML", "B_1_filter.idXML"]


def test_every_file_kept_when_percolator_succeeds(tmp_path, monkeypatch):
    wf, names = _workflow(tmp_path, monkeypatch, percolator_fails=[])
    wf.execution()
    lfq = _inputs(wf, "ProteomicsLFQ")
    assert [Path(f).name for f in lfq["in"][0]] == names


def test_run_stops_when_percolator_fails_on_every_file(tmp_path, monkeypatch):
    wf, _ = _workflow(tmp_path, monkeypatch, percolator_fails=[0, 1, 2])
    assert wf.execution() is False
    assert not any(t == "IDFilter" for t, _, _ in wf.executor.calls)


def test_rescoring_all_files_together_runs_one_command(tmp_path, monkeypatch):
    """With 'Rescore all files together' one PercolatorAdapter call takes every
    search result; its output is split back into one result per file."""
    import src.common.pooled_fdr as pooled_fdr

    wf, names = _workflow(tmp_path, monkeypatch, percolator_fails=[])
    wf.params["rescore-pooled"] = True
    split = []
    monkeypatch.setattr(
        pooled_fdr, "split_pooled",
        lambda merged, ins, outs: split.append((merged, list(ins), list(outs))) or [1] * len(ins),
    )
    wf.execution()

    percolator = [c for c in wf.executor.calls if c[0] == "PercolatorAdapter"]
    assert len(percolator) == 1
    io = percolator[0][1]
    assert len(io["in"]) == 1 and len(io["in"][0]) == 3  # one command, three inputs
    assert io["out"] == [str(tmp_path / "results" / "percolator_all_files" / "all_files_per.idXML")]

    merged, ins, outs = split[0]
    assert merged == io["out"][0] and ins == io["in"][0]
    assert [Path(f).name for f in outs] == ["A_1_per.idXML", "Blank_per.idXML", "B_1_per.idXML"]
    idfilter = _inputs(wf, "IDFilter")
    assert [Path(f).name for f in idfilter["in"]] == [Path(f).name for f in outs]


def test_one_file_is_rescored_on_its_own(tmp_path, monkeypatch):
    wf, _ = _workflow(tmp_path, monkeypatch, percolator_fails=[])
    wf.params["rescore-pooled"] = True
    assert wf.rescores_pooled(["only.idXML"]) is False
    assert wf.rescores_pooled(["a.idXML", "b.idXML"]) is True
