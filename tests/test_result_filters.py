"""Filters of the PSM and protein tables."""
import pandas as pd

from src.common.result_filters import (
    add_protein_flags, file_conditions, max_replicates, protein_numeric_filters,
    q_choices, q_suggestions, replicates_found,
)


def test_q_choices_never_loosen_the_run():
    assert q_choices(1.0) == [1.0, 0.5, 0.1, 0.01]
    assert q_choices(5.0) == [5.0, 1.0, 0.5, 0.1, 0.01]
    assert q_choices(0.1) == [0.1, 0.01]
    assert q_choices(2.0) == [2.0, 1.0, 0.5, 0.1, 0.01]
    # PSM filter off: any level, the first meaning no threshold
    assert q_choices(None) == [100.0, 5.0, 1.0, 0.5, 0.1, 0.01]


def test_q_suggestions_start_at_the_run():
    assert q_suggestions(1.0) == [
        {"value": 0.01, "label": "1 % (run)"},
        {"value": 0.005, "label": "0.5 %"},
        {"value": 0.001, "label": "0.1 %"},
        {"value": 0.0001, "label": "0.01 %"},
    ]
    # PSM filter off: the standard levels; "no threshold" is the dialog's "Any"
    assert [s["value"] for s in q_suggestions(None)] == [0.05, 0.01, 0.005, 0.001, 0.0001]


def test_protein_numeric_filters():
    assert protein_numeric_filters(0) == {"peptides": {"op": ">=", "suggestions": [2, 3, 5]}}
    assert protein_numeric_filters(3)["replicates"] == {"op": ">=", "suggestions": [1, 2, 3]}


def _long():
    # P1 in both files of A and one of B; P2 in one file of each; CON__ in one
    return pd.DataFrame({
        "accession": ["P1", "P1", "P1", "P1", "P2", "P2", "CON__sp|K1C10_HUMAN|"],
        "file": ["a1", "a1", "a2", "b1", "a1", "b2", "b1"],
    })


def test_replicates_found_is_the_best_condition():
    conditions = {"a1": "A", "a2": "A", "b1": "B", "b2": "B", "x": ""}
    counts = replicates_found(_long(), conditions)
    assert counts["P1"] == 2 and counts["P2"] == 1
    assert max_replicates(file_conditions(["a1", "a2", "b1", "c"], lambda f: {"a1": "A", "a2": "A", "b1": "B"}.get(f, ""))) == 2
    assert replicates_found(_long(), {}).empty


def test_protein_flags():
    proteins = pd.DataFrame({
        "accession": ["P1", "P2", "CON__sp|K1C10_HUMAN|"],
        "gene": ["ALB", "KRT1", ""],
        "found in": ["3 of 4", "2 of 4", "1 of 4"],
        "description": ["Serum albumin", "Keratin", "Keratin, type I"],
        "peptides": [5, 1, 2],
    })
    flagged = add_protein_flags(proteins, _long(), {"a1": "A", "a2": "A", "b1": "B", "b2": "B"})
    assert flagged["contaminant"].tolist() == [False, False, True]
    assert flagged["replicates"].tolist() == [2, 1, 1]
    assert list(flagged.columns).index("replicates") == list(flagged.columns).index("found in") + 1
    # without conditions there is no replicates column
    assert "replicates" not in add_protein_flags(proteins, _long(), {}).columns


def test_replicates_count_samples_not_fractions():
    from src.common.result_filters import max_replicates, replicates_found

    psms = pd.DataFrame({"accession": ["P1", "P1", "P1"], "file": ["s1_f1.mzML", "s1_f2.mzML", "s2_f1.mzML"]})
    conditions = {"s1_f1.mzML": "A", "s1_f2.mzML": "A", "s2_f1.mzML": "A", "s2_f2.mzML": "A"}
    sample_of = {"s1_f1.mzML": "A_1", "s1_f2.mzML": "A_1", "s2_f1.mzML": "A_2", "s2_f2.mzML": "A_2"}
    # per file, three fractions look like three replicates; per sample it is two
    assert replicates_found(psms, conditions)["P1"] == 3
    assert replicates_found(psms, conditions, sample_of)["P1"] == 2
    assert max_replicates(conditions) == 4
    assert max_replicates(conditions, sample_of) == 2
