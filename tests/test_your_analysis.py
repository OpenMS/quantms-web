"""Your analysis on Describe experiment: what a start does, and vendor runs keeping their groups."""
import json
from types import SimpleNamespace

from src.WorkflowTest import CLOSED_SEARCH, OPEN_SEARCH, WorkflowTest


def test_summary_names_the_analysis():
    text = WorkflowTest.summary("LFQ", 6, "Orbitrap HCD", CLOSED_SEARCH, "yeast.fasta")
    assert text.startswith("**A label-free analysis of 6 runs**")
    assert "yeast.fasta" in text and "Orbitrap (HCD) tolerances" in text
    text = WorkflowTest.summary("TMT", 1, None, OPEN_SEARCH, None)
    assert text.startswith("**An isobaric (TMT/iTRAQ) analysis of 1 run**")
    assert "an open database search" in text and "your own tolerances" in text
    assert "once you choose a protein database" in text


def test_summary_without_runs_asks_for_them():
    assert "Browse" in WorkflowTest.summary("LFQ", 0, None, CLOSED_SEARCH, None)


def test_converted_runs_keep_their_groups(tmp_path):
    params_file = tmp_path / "params.json"
    params = {"mzML-files": ["/in/a.raw", "/in/b.mzML"], "mzML-group-a.raw": "ctrl", "mzML-group-b.mzML": "drug"}
    params_file.write_text(json.dumps(params))
    wf = object.__new__(WorkflowTest)
    wf.params = dict(params)
    wf.parameter_manager = SimpleNamespace(
        params_file=params_file,
        get_parameters_from_json=lambda: json.loads(params_file.read_text()),
    )
    wf.alias_converted_runs(["/in/a.raw", "/in/b.mzML"], ["/wf/converted-mzML/a.mzML", "/in/b.mzML"])
    stored = json.loads(params_file.read_text())
    assert stored["mzML-group-a.mzML"] == "ctrl"
    assert stored["mzML-group-b.mzML"] == "drug"
    assert wf.params["mzML-group-a.mzML"] == "ctrl"


def test_summary_names_the_labels_and_the_library():
    text = WorkflowTest.summary("TMT", 4, "Orbitrap HCD", CLOSED_SEARCH, "human.fasta", plex="tmt16plex")
    assert text.startswith("**An isobaric (TMT/iTRAQ) analysis of 4 runs** labelled with TMT 16")
    text = WorkflowTest.summary("LFQ", 2, None, CLOSED_SEARCH, None, library=True)
    assert text.endswith("and a spectral library of the identified peptides.")


def test_steps_follow_the_analysis():
    wf = object.__new__(WorkflowTest)
    wf.mode = lambda: "LFQ"
    wf.library = lambda: False
    wf.ptm_enrichment_on = lambda: False
    assert WorkflowTest.steps(wf) == ["Samples and groups", "Database search", "Rescoring",
                                      "False discovery rate", "Quantification"]
    wf.library = lambda: True
    assert WorkflowTest.steps(wf)[-1] == "Spectral library"
    wf.mode = lambda: "TMT"
    assert WorkflowTest.steps(wf)[1] == "Labels and decoys"
    wf.ptm_enrichment_on = lambda: True
    assert WorkflowTest.steps(wf)[-1] == "PTM enrichment"


def test_summary_names_the_ptm_enrichment():
    text = WorkflowTest.summary("LFQ", 4, None, CLOSED_SEARCH, None, ptm=True)
    assert text.endswith("and the modifications enriched among the changed proteins.")
    text = WorkflowTest.summary("LFQ", 4, None, CLOSED_SEARCH, None, library=True, ptm=True)
    assert "modifications enriched among the changed proteins and a spectral library" in text
