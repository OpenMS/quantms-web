"""A file the search fails on is left out; the rest of the study goes on.

ProSE searches several files in one process, so one file it cannot read
failed the whole batch and the run. run_search() now asks the executor to run
a failed batch again file by file (``retry_singly``, where the executor has
it), and execution() drops the files that still fail from every list the next
steps take.
"""
from types import SimpleNamespace

# Imported inside the tests, as in test_search_batches.py.


class RetryExecutor:
    """An executor with ``retry_singly``: the second file fails."""

    def __init__(self, fails=(1,)):
        self.fails = list(fails)
        self.failed = []
        self.calls = []

    def _get_max_threads(self):
        return 16

    def run_topp(self, tool, io, params=None, tool_instance_name=None,
                 allow_failed=False, retry_singly=False):
        self.calls.append((tool, io, retry_singly))
        self.failed = list(self.fails) if retry_singly else []
        return True

    def keep_succeeded(self, values):
        return [v for i, v in enumerate(values) if i not in set(self.failed)]


class OldExecutor:
    """An executor from before ``retry_singly``."""

    def __init__(self):
        self.failed = [0]  # left over from an earlier step
        self.calls = []

    def _get_max_threads(self):
        return 16

    def run_topp(self, tool, io, params=None, tool_instance_name=None, allow_failed=False):
        self.calls.append((tool, io))
        return True


def _wf(executor):
    from src.WorkflowTest import WorkflowTest

    wf = object.__new__(WorkflowTest)
    wf.executor = executor
    wf.logger = SimpleNamespace(log=lambda *a, **k: None)
    return wf


def test_run_search_asks_for_a_retry_file_by_file(tmp_path):
    wf = _wf(RetryExecutor())
    assert wf.run_search("ProSE", ["a.mzML", "b.mzML"], ["a.idXML", "b.idXML"], tmp_path / "db.fasta")
    tool, io, retry = wf.executor.calls[0]
    assert tool == "ProSE" and retry is True
    assert io["in"] == [["a.mzML", "b.mzML"]]
    assert wf.executor.failed == [1]


def test_run_search_with_an_older_executor_keeps_every_file(tmp_path):
    wf = _wf(OldExecutor())
    assert wf.run_search("ProSE", ["a.mzML", "b.mzML"], ["a.idXML", "b.idXML"], tmp_path / "db.fasta")
    assert wf.executor.calls[0][0] == "ProSE"
    assert wf.executor.failed == []  # nothing from an earlier step is dropped


def test_supports_retry_singly():
    from src.WorkflowTest import supports_retry_singly

    assert supports_retry_singly(RetryExecutor())
    assert not supports_retry_singly(OldExecutor())
    assert not supports_retry_singly(SimpleNamespace(run_topp=None))


def test_the_file_the_search_fails_on_is_left_out(tmp_path, monkeypatch):
    """Percolator gets only the files the search read (LFQ)."""
    from tests import test_percolator_left_out as perc

    wf, names = perc._workflow(tmp_path, monkeypatch, percolator_fails=[])

    def search(instance, in_mzml, out_idxml, fasta):
        wf.executor.failed = [1]  # Blank.mzML
        return True

    monkeypatch.setattr(wf, "run_search", search, raising=False)
    assert wf.execution() is False  # the fake ProteomicsLFQ stops the run

    rescored = perc._inputs(wf, "PercolatorAdapter")["in"]
    assert [p.split("/")[-1].split("_prose")[0].split(".")[0] for p in rescored] == ["A_1", "B_1"]
