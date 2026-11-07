"""Views a run cannot fill: the group comparison needs samples in two groups."""
import json
from pathlib import Path

import streamlit as st

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


def _needs(component):
    need = component.get("needs") or []
    return {need} if isinstance(need, str) else set(need)


def test_components_comparing_samples_need_samples_or_groups():
    needs = {c["name"]: _needs(c) for c in CATALOG}
    assert {n for n, v in needs.items() if analyses.GROUPS in v} == {
        "group_comparison", "volcano", "go_enrichment", "comparison_missing_values",
        "comparison_imputed_values", "comparison_normalization", "comparison_pvalues",
    }
    assert {n for n, v in needs.items() if analyses.SAMPLES in v} == {"heatmap", "clustered_heatmap"}
    assert set().union(*needs.values()) == set(analyses.ANALYSES)


def test_every_component_needs_the_results_it_draws():
    needs = {c["name"]: _needs(c) for c in CATALOG}
    assert {n for n, v in needs.items() if analyses.QUANTIFICATION in v} == {
        "protein_abundance", "pca", "heatmap", "clustered_heatmap", "group_comparison", "volcano", "go_enrichment",
        "comparison_missing_values", "comparison_imputed_values", "comparison_normalization", "comparison_pvalues"}
    assert needs["psms_rescored"] == {analyses.RESCORING}
    assert needs["psms_search"] == needs["id_quality"] == {analyses.SEARCH}
    assert all(needs[n] for n in needs)


def _results(workspace, run_id, folder, name):
    from src.common.results_helpers import get_results_dir
    directory = get_results_dir(workspace, run_id) / folder
    directory.mkdir(parents=True, exist_ok=True)
    (directory / name).write_text("x")


def test_results_lock_what_the_picked_files_lack(tmp_path):
    first = _run(tmp_path, "First", {"mzML-files": ["/in/a.mzML"]})
    assert analyses.available(tmp_path, first) == set()
    _results(tmp_path, first, "psm_filter", "a.idXML")
    _results(tmp_path, first, "search_results", "a.idXML")
    assert analyses.available(tmp_path, first) == {analyses.IDENTIFICATION, analyses.SEARCH}

    # Abundances of any run of the experiment count; identifications those of the picked run
    added = _added_run(tmp_path, first, "Added", {"mzML-files": ["/in/b.mzML"]})
    _results(tmp_path, added, "quant_results", "q.csv")
    assert analyses.QUANTIFICATION in analyses.available(tmp_path, first)
    st.session_state["viewer-member-" + first + "-e0"] = added
    try:
        has = analyses.available(tmp_path, first, 0)
    finally:
        st.session_state.pop("viewer-member-" + first + "-e0")
    assert analyses.IDENTIFICATION not in has and analyses.QUANTIFICATION in has


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


def test_options_only_process_more_files(tmp_path):
    single = _run(tmp_path, "One file", {"mzML-files": ["/in/a.mzML"]})
    two = _run(tmp_path, "Two, no groups", {"mzML-files": ["/in/d.mzML", "/in/e.mzML"]})
    for run_id in (single, two):
        labels = [o["label"] for o in analyses.options(tmp_path, run_id, [analyses.GROUPS])]
        assert labels == ["Process more files with these settings"]


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

    # The experiment is listed once
    runs = RunManager(Path(tmp_path, "topp-workflow"))
    assert [r["id"] for r in runs.experiments()] == [first]
