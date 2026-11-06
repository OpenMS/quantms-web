"""Views a run cannot fill: the group comparison needs samples in two groups."""
import json
from pathlib import Path

from src.viewer import analyses
from src.viewer.catalog import CATALOG
from src.workflow.RunManager import FINISHED, RunManager


def _run(workspace, name, params):
    """A finished run started with ``params`` (the workflow's params.json then)."""
    workflow_dir = Path(workspace, "topp-workflow")
    workflow_dir.mkdir(parents=True, exist_ok=True)
    params_file = workflow_dir / "params.json"
    params_file.write_text(json.dumps(params))
    runs = RunManager(workflow_dir)
    run_id = runs.new_run(name, params_file)
    runs.finish(run_id, FINISHED)
    return run_id


def _set_current(workspace, params):
    Path(workspace, "topp-workflow", "params.json").write_text(json.dumps(params))


def test_components_comparing_samples_need_samples_or_groups():
    needs = {c["name"]: c.get("needs") for c in CATALOG}
    assert {n for n, v in needs.items() if v == analyses.GROUPS} == {"group_comparison", "volcano", "go_enrichment"}
    assert {n for n, v in needs.items() if v == analyses.SAMPLES} == {"heatmap", "clustered_heatmap"}
    assert set(needs.values()) - {None} == set(analyses.ANALYSES)


def test_a_run_has_groups_once_two_of_its_samples_differ(tmp_path):
    files = ["/in/ctrl.mzML", "/in/treated.mzML"]
    run_id = _run(tmp_path, "Two files", {"mzML-files": files})
    assert analyses.available(tmp_path, run_id) == {analyses.SAMPLES}
    assert analyses.run_samples(tmp_path, run_id) == 2
    one = _run(tmp_path, "One file", {"mzML-files": ["/in/a.mzML"]})
    assert analyses.available(tmp_path, one) == set()

    # Groups set on Describe experiment after the run apply to it
    _set_current(tmp_path, {"mzML-group-ctrl.mzML": "control", "mzML-group-treated.mzML": "control"})
    assert analyses.available(tmp_path, run_id) == {analyses.SAMPLES}
    _set_current(tmp_path, {"mzML-group-ctrl.mzML": "control", "mzML-group-treated.mzML": "treated"})
    assert analyses.available(tmp_path, run_id) == {analyses.SAMPLES, analyses.GROUPS}
    assert analyses.available(tmp_path, None) is None


def test_tmt_channels_skipped_take_no_part(tmp_path):
    run_id = _run(tmp_path, "TMT", {"analysis-mode": "TMT"})
    _set_current(tmp_path, {"TMT-group-sample1": "A", "TMT-group-sample2": "skip"})
    assert analyses.GROUPS not in analyses.available(tmp_path, run_id)
    _set_current(tmp_path, {"TMT-group-sample1": "A", "TMT-group-sample2": "B"})
    assert analyses.GROUPS in analyses.available(tmp_path, run_id)


def test_options_add_groups_or_files_and_bring_in_a_run_with_groups(tmp_path):
    single = _run(tmp_path, "One file", {"mzML-files": ["/in/a.mzML"]})
    grouped = _run(tmp_path, "Grouped", {"mzML-files": ["/in/b.mzML", "/in/c.mzML"]})
    _set_current(tmp_path, {"mzML-group-b.mzML": "x", "mzML-group-c.mzML": "y"})
    args = dict(selection_key=lambda i: f"sel{i}", view_args=None)

    labels = [o["label"] for o in analyses.options(tmp_path, single, [analyses.GROUPS], None, **args)]
    # One file has nothing to group: more files are the way
    assert labels == ["Process more files with these settings", "Show Grouped",
                      "Compare side by side with Grouped"]

    two = _run(tmp_path, "Two, no groups", {"mzML-files": ["/in/d.mzML", "/in/e.mzML"]})
    labels = [o["label"] for o in analyses.options(tmp_path, two, [analyses.GROUPS], None, index=1, **args)]
    assert labels[:2] == ["Give the samples groups", "Process more files with these settings"]
    assert not any(label.startswith("Compare") for label in labels)
    assert grouped != two
