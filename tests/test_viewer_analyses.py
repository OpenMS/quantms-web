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


def _added_run(workspace, first, name, params):
    """A finished run that processed more files for ``first``'s experiment."""
    workflow_dir = Path(workspace, "topp-workflow")
    params_file = workflow_dir / "params.json"
    params_file.write_text(json.dumps(params))
    runs = RunManager(workflow_dir)
    run_id = runs.new_run(name, params_file, adds_to=first)
    runs.finish(run_id, FINISHED)
    return run_id


def _quant(workspace, run_id, rows):
    quant = RunManager(Path(workspace, "topp-workflow")).results_dir(run_id) / "quant_results"
    quant.mkdir(parents=True, exist_ok=True)
    lines = ["ProteinName,PeptideSequence,Reference,Intensity"]
    lines += [f"{p},{s},{r},{i}" for p, s, r, i in rows]
    (quant / "msstats.csv").write_text("\n".join(lines) + "\n")


def test_more_files_join_the_experiment_and_its_comparison(tmp_path):
    from src.common.results_helpers import load_abundance_data

    first = _run(tmp_path, "Controls", {"mzML-files": ["/in/ctrl1.mzML", "/in/ctrl2.mzML"],
                                         "mzML-group-ctrl1.mzML": "control",
                                         "mzML-group-ctrl2.mzML": "control"})
    # All of one group: nothing to compare yet
    assert analyses.available(tmp_path, first) == {analyses.SAMPLES}
    added = _added_run(tmp_path, first, "Treated", {"mzML-files": ["/in/drug1.mzML"],
                                                      "mzML-group-drug1.mzML": "treated"})
    # Configure has since moved on to other files: each run keeps the groups it set
    _set_current(tmp_path, {})
    assert analyses.available(tmp_path, first) == {analyses.SAMPLES, analyses.GROUPS}
    assert analyses.run_samples(tmp_path, first) == 3
    assert analyses.run_groups(tmp_path, first) == {
        "ctrl1.mzML": "control", "ctrl2.mzML": "control", "drug1.mzML": "treated"}

    _quant(tmp_path, first, [("P1", "PEPTIDE", "ctrl1.mzML", 10), ("P1", "PEPTIDE", "ctrl2.mzML", 12),
                             ("P2", "OTHER", "ctrl1.mzML", 5)])
    _quant(tmp_path, added, [("P1", "PEPTIDE", "drug1.mzML", 40), ("P3", "NEW", "drug1.mzML", 7)])
    pivot, _, group_map = load_abundance_data(str(tmp_path), 0.0, 0.0, first, (first, added))
    assert list(pivot.columns) == ["ProteinName", "ctrl1", "ctrl2", "drug1", "PeptideSequence"]
    p1 = pivot.set_index("ProteinName").loc["P1"]
    assert (p1["ctrl1"], p1["ctrl2"], p1["drug1"]) == (10, 12, 40)
    assert set(pivot["ProteinName"]) == {"P1", "P2", "P3"}
    assert group_map["drug1.mzML"] == "treated"

    # The experiment is listed once, and "Show" offers no run of it
    runs = RunManager(Path(tmp_path, "topp-workflow"))
    assert [r["id"] for r in runs.experiments()] == [first]
    args = dict(selection_key=lambda i: f"sel{i}", view_args=None)
    labels = [o["label"] for o in analyses.options(tmp_path, first, [analyses.GROUPS], None, **args)]
    assert not any(label.startswith("Show") for label in labels)
