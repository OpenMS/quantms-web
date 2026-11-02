"""
A local workflow run must not inherit the server's thread pools.

The server process uses polars (and pyOpenMS's OpenMP pool) on the results
pages. A workflow process forked from it got those pools without their
threads, and the run blocked forever at its first polars query: in the desktop
app every run started after a results page was open never finished.

This test warms polars' pool in the test process, the way a results page does,
then starts a workflow whose execution() runs a polars query, and checks that
the run finishes.
"""

import sys
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

pl = pytest.importorskip("polars")
pytest.importorskip("streamlit")
pytest.importorskip("pyopenms")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.workflow.WorkflowManager import WorkflowManager  # noqa: E402


class PolarsWorkflow(WorkflowManager):
    def execution(self) -> bool:
        df = pl.DataFrame({"g": [i % 7 for i in range(200_000)], "v": range(200_000)})
        groups = df.lazy().group_by("g").agg(pl.col("v").sum()).collect()
        self.logger.log(f"groups: {groups.height}")
        return True


def _warm_polars_pool():
    df = pl.DataFrame({"g": [i % 100 for i in range(2_000_000)], "v": range(2_000_000)})
    df.lazy().group_by("g").agg(pl.col("v").sum()).sort("g").collect()


def test_run_finishes_after_the_server_used_polars(tmp_path, monkeypatch):
    monkeypatch.chdir(ROOT)  # the workflow process imports src.* from the app folder
    _warm_polars_pool()

    workflow_dir = tmp_path / "polars-workflow"
    workflow_dir.mkdir()
    wf = object.__new__(PolarsWorkflow)  # skips __init__, which needs a Streamlit session
    wf.workflow_dir = workflow_dir
    wf.run_id = None
    wf.executor = SimpleNamespace(pid_dir=workflow_dir / "pids")
    # The test module is imported in the workflow process by this name
    monkeypatch.setattr(PolarsWorkflow, "__module__", "tests.test_workflow_local_process")

    wf._start_workflow_local()

    log = workflow_dir / "logs" / "minimal.log"
    deadline = time.time() + 120
    while time.time() < deadline:
        if log.exists() and "WORKFLOW FINISHED" in log.read_text(encoding="utf-8"):
            break
        time.sleep(0.5)
    text = log.read_text(encoding="utf-8") if log.exists() else ""
    assert "groups: 7" in text and "WORKFLOW FINISHED" in text, f"run did not finish:\n{text}"
    for _ in range(20):  # the pid dir goes once the process has cleaned up
        if not (workflow_dir / "pids").exists():
            break
        time.sleep(0.5)
    assert not (workflow_dir / "pids").exists()
