"""Tests for src/common/fractions.py and the results that merge a sample's fractions."""
# Sorted after test_parameter_defaults on purpose: the samples and protein-table
# tests import modules that pull in src.workflow.ParameterManager.
import pandas as pd

from src.common.fractions import design_samples, fraction_samples, sample_of
from src.workflow.ExperimentalDesign import save_design


def _design(rows):
    return pd.DataFrame(rows, columns=["file", "condition", "replicate", "fraction"])


FRACTIONATED = _design([
    ("a1.mzML", "ctrl", 1, 1), ("a2.mzML", "ctrl", 1, 2),
    ("b1.mzML", "ctrl", 2, 1), ("b2.mzML", "ctrl", 2, 2),
    ("c1.mzML", "treated", 1, 1), ("c2.mzML", "treated", 1, 2),
])


def test_design_samples_names_samples_by_condition_and_replicate():
    samples = design_samples(FRACTIONATED, {})
    assert {f: s["sample"] for f, s in samples.items()} == {
        "a1.mzML": "ctrl_1", "a2.mzML": "ctrl_1",
        "b1.mzML": "ctrl_2", "b2.mzML": "ctrl_2",
        "c1.mzML": "treated_1", "c2.mzML": "treated_1",
    }


def test_design_samples_takes_conditions_from_the_run_parameters():
    samples = design_samples(FRACTIONATED, {"mzML-group-c1.mzML": "drug"})
    assert samples["c2.mzML"] == {"sample": "drug_1", "condition": "drug"}


def test_technical_replicates_stay_apart():
    design = _design([
        ("t1f1.mzML", "ctrl", 1, 1), ("t1f2.mzML", "ctrl", 1, 2),
        ("t2f1.mzML", "ctrl", 1, 1), ("t2f2.mzML", "ctrl", 1, 2),
    ])
    names = {f: s["sample"] for f, s in design_samples(design, {}).items()}
    assert names == {"t1f1.mzML": "ctrl_1_t1", "t1f2.mzML": "ctrl_1_t1", "t2f1.mzML": "ctrl_1_t2", "t2f2.mzML": "ctrl_1_t2"}


def test_fraction_samples_only_for_fractionated_runs(tmp_path):
    assert fraction_samples(tmp_path, {}) is None  # no design passed to ProteomicsLFQ
    quant = tmp_path / "quant_results"
    quant.mkdir()
    save_design(quant, FRACTIONATED.assign(fraction=1))
    assert fraction_samples(tmp_path, {}) is None
    save_design(quant, FRACTIONATED)
    assert fraction_samples(tmp_path, {})["b2.mzML"]["sample"] == "ctrl_2"


def test_sample_of_falls_back_to_the_file_stem():
    samples = design_samples(FRACTIONATED, {})
    assert sample_of("a2.mzML", samples) == "ctrl_1"
    assert sample_of("other.mzML", samples) == "other"
    assert sample_of("x.mzML", None) == "x"


def test_lfq_samples_merge_fractions(tmp_path):
    from src.common.samples import lfq_samples
    from tests.test_samples import _idxml, _quant_csv

    results = tmp_path / "results"
    _idxml(results / "psm_filter" / "a1_filter.idXML", [("PEPTIDE", "P1", "target")])
    _idxml(results / "psm_filter" / "a2_filter.idXML", [("PEPTIDE", "P1", "target"), ("PEPTIDEK", "P2", "target")])
    _idxml(results / "psm_filter" / "c1_filter.idXML", [("PEPTIDE", "P1", "target")])
    _quant_csv(results / "quant_results" / "q.csv", [
        ("P1", "PEPTIDE", "a1.mzML", 1.0), ("P2", "PEPTIDEK", "a2.mzML", 1.0), ("P1", "PEPTIDE", "c1.mzML", 1.0),
    ])
    save_design(results / "quant_results", _design([
        ("a1.mzML", "ctrl", 1, 1), ("a2.mzML", "ctrl", 1, 2), ("c1.mzML", "treated", 1, 1), ("c2.mzML", "treated", 1, 2),
    ]))
    groups = {"a1.mzML": "ctrl", "a2.mzML": "ctrl", "c1.mzML": "treated", "c2.mzML": "treated"}
    df = lfq_samples(tmp_path, ["in/a1.mzML", "in/a2.mzML", "in/c1.mzML", "in/c2.mzML"], groups)
    assert df["Sample"].tolist() == ["ctrl_1", "treated_1"]
    ctrl = df.iloc[0]
    # PSMs add up over fractions, peptides and proteins are counted once
    assert (ctrl["Target PSMs"], ctrl["Peptides"], ctrl["Proteins"]) == (3, 2, 2)
    assert (ctrl["Quantified proteins"], ctrl["Quantified peptides"]) == (2, 2)
    assert df.iloc[1]["Target PSMs"] == 1


def test_protein_table_counts_per_sample():
    from src.common.protein_table import build_tables

    psms = pd.DataFrame({
        "file": ["a1.mzML", "a2.mzML", "c1.mzML"],
        "scan": [1, 2, 3], "rt": [1.0, 2.0, 3.0], "mz": [400.0, 400.0, 400.0], "charge": [2, 2, 2],
        "peptide": ["PEPTIDEK"] * 3, "sequence": ["PEPTIDEK"] * 3,
        "score": [0.1] * 3, "q_value": [0.001] * 3, "proteins": [["P1"]] * 3,
    })
    samples = {"a1.mzML": "ctrl_1", "a2.mzML": "ctrl_1", "c1.mzML": "treated_1"}
    proteins, peptides, _, _ = build_tables(psms, {}, ["a1.mzML", "a2.mzML", "c1.mzML"], samples)
    p1 = proteins.iloc[0]
    assert p1["found in"] == "2 of 2"
    assert (p1["PSMs ctrl_1"], p1["PSMs treated_1"]) == (2, 1)
    assert peptides.iloc[0]["found in samples"] == 2
