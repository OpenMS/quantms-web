"""Filter bars on the PSM and protein pages."""
import pandas as pd

from src.common.result_filters import (
    add_protein_flags, file_conditions, filter_proteins, filter_psms_by_q, max_replicates,
    q_choices, replicates_found, text_mask,
)


def test_q_choices_never_loosen_the_run():
    assert q_choices(1.0) == [1.0, 0.5, 0.1, 0.01]
    assert q_choices(5.0) == [5.0, 1.0, 0.5, 0.1, 0.01]
    assert q_choices(0.1) == [0.1, 0.01]
    assert q_choices(2.0) == [2.0, 1.0, 0.5, 0.1, 0.01]
    # PSM filter off: any level, the first meaning no threshold
    assert q_choices(None) == [100.0, 5.0, 1.0, 0.5, 0.1, 0.01]


def test_text_mask_matches_every_word_in_any_column():
    df = pd.DataFrame({"a": ["ALBU_HUMAN", "KRT1", None], "b": ["Serum albumin", "Keratin", "x"]})
    assert text_mask(df, ["a", "b"], "albu").tolist() == [True, False, False]
    assert text_mask(df, ["a", "b"], "serum ALBU").tolist() == [True, False, False]
    assert text_mask(df, ["a", "b"], "serum keratin").tolist() == [False, False, False]
    assert text_mask(df, ["a", "b"], "").tolist() == [True, True, True]
    # literal, not a regular expression
    assert text_mask(df, ["a"], "KRT1(").tolist() == [False, False, False]


def test_filter_psms_by_q():
    psms = pd.DataFrame({"q_value": [0.0001, 0.004, 0.009, float("nan")]})
    assert len(filter_psms_by_q(psms, None)) == 4
    assert filter_psms_by_q(psms, 0.005)["q_value"].tolist()[:2] == [0.0001, 0.004]
    assert len(filter_psms_by_q(psms, 0.005)) == 3  # no q-value: kept


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


def test_protein_filters():
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
    assert filter_proteins(flagged, "keratin")["accession"].tolist() == ["P2", "CON__sp|K1C10_HUMAN|"]
    assert filter_proteins(flagged, min_peptides=2)["accession"].tolist() == ["P1", "CON__sp|K1C10_HUMAN|"]
    assert filter_proteins(flagged, min_replicates=2)["accession"].tolist() == ["P1"]
    # without conditions there is no replicates column, and the threshold is ignored
    plain = add_protein_flags(proteins, _long(), {})
    assert "replicates" not in plain.columns
    assert len(filter_proteins(plain, min_replicates=2)) == 3
