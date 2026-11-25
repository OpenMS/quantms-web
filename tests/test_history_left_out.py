"""History of a result: the files a run went on without."""
import sys
import types

from src.viewer import history_ui


class _Runs:
    """A RunManager with left-out files (the template's left_out_files)."""

    def __init__(self, members, left_out):
        self._members = members
        self._left_out = left_out

    def experiment(self, run_id):
        return self._members

    def experiment_root(self, run_id):
        return self._members[0]

    def left_out_files(self, run_id, key):
        assert key == "mzML-files"
        return self._left_out.get(run_id, [])


def test_left_out_files_cover_every_run_of_the_experiment():
    runs = _Runs(["a", "b"], {"a": ["/d/x.mzML"], "b": ["/d/y.mzML", "/d/x.mzML"]})
    assert history_ui.left_out_files(runs, "b") == ["/d/x.mzML", "/d/y.mzML"]


def test_left_out_files_without_the_template_feature():
    class Old:
        pass

    assert history_ui.left_out_files(Old(), "a") == []


def test_process_left_out_adds_the_files_to_the_experiment(monkeypatch):
    calls = []

    class FakeWorkflow:
        def follow_up(self, run_id, **kwargs):
            calls.append((run_id, kwargs))
            return True

    monkeypatch.setitem(sys.modules, "src.WorkflowTest", types.SimpleNamespace(WorkflowTest=FakeWorkflow))
    history_ui.process_left_out(_Runs(["a", "b"], {}), "b", ["/d/x.mzML"])
    run_id, kwargs = calls[0]
    assert run_id == "a"
    assert kwargs["add_files"] == "mzML-files" and kwargs["files"] == ["/d/x.mzML"]
    assert kwargs["page"] == history_ui.CONFIGURE_PAGE
