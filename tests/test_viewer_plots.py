"""Tables behind the Viewer's Insight plots, GO enrichment of picked proteins and the change cut-offs."""
import numpy as np
import pandas as pd

from src.common.differential import CHANGE_DEFAULTS, change_cutoffs
from src.common.postprocessing_plots import (
    imputed_value_preview, missingness_vs_intensity, pvalue_histogram, sample_distributions,
)
from src.viewer import go_terms


def test_change_cutoffs_default_and_saved():
    assert change_cutoffs({}) == (CHANGE_DEFAULTS["postproc-change-fc"], CHANGE_DEFAULTS["postproc-change-padj"])
    assert change_cutoffs({"postproc-change-fc": 0.5, "postproc-change-padj": "0.01"}) == (0.5, 0.01)
    assert change_cutoffs({"postproc-change-fc": "x"}) == (1.0, 0.05)


def test_pvalue_histogram_bins_and_expected_count():
    table, expected = pvalue_histogram(pd.DataFrame({"p-value": [0.01, 0.02, 0.5, 0.99, None]}))
    assert len(table) == 20 and table["Proteins"].sum() == 4
    assert table["Proteins"].iloc[0] == 2 and table["Proteins"].iloc[-1] == 1
    assert expected == 4 / 20


def test_missing_and_imputed_tables():
    before = pd.DataFrame({"s1": [8.0, 16.0, 0.0], "s2": [8.0, 16.0, 4.0]})
    after = pd.DataFrame({"s1": [8.0, 16.0, 2.0], "s2": [8.0, 16.0, 4.0]})
    table, gaps = missingness_vs_intensity(before, ["s1", "s2"])
    assert gaps == 1 and set(table["series"]) == {"complete", "missing"}
    filled = imputed_value_preview(before, after, ["s1", "s2"])
    assert filled.loc[filled["series"] == "imputed", "y"].sum() == 1
    assert filled.loc[filled["series"] == "observed", "y"].sum() == 5


def test_sample_distributions_keep_quantiles():
    rng = np.random.default_rng(0)
    before = pd.DataFrame({"s1": 2 ** rng.normal(20, 1, 1000)})
    after = np.log2(before)
    table = sample_distributions(before, after, ["s1"], after_is_log=True, points=101)
    assert set(table["Stage"]) == {"before", "after"}
    kept = table[table["Stage"] == "after"]["Intensity"]
    assert len(kept) == 101
    assert abs(kept.median() - after["s1"].median()) < 1e-9


def test_go_terms_of_picked_proteins(tmp_path):
    proteins = ["sp|P1|A_HUMAN", "sp|P2|B_HUMAN", "sp|P3|C_HUMAN", "sp|P4|D_HUMAN"]
    calls = []

    def fetch(accessions):
        calls.append(accessions)
        return pd.DataFrame({
            "query": ["P1", "P2", "P3", "P4"],
            "BP_terms": [["kinase activity"], ["kinase activity"], ["transport"], None],
            "CC_terms": [None] * 4,
            "MF_terms": [None] * 4,
        })

    index, error = go_terms.term_index(tmp_path, proteins, fetch=fetch)
    assert error is None and go_terms.index_ready(tmp_path, proteins)
    go_terms.term_index(tmp_path, proteins, fetch=fetch)
    assert len(calls) == 1  # looked up once per protein set
    annotated, counts = go_terms.term_counts(index, proteins)
    assert annotated == 3
    assert counts["BP"]["term"].to_list()[0] == "kinase activity"
    result = go_terms.enrichment(index, proteins[:2], proteins)
    table = go_terms.term_table(result["categories"]["BP"], {"P1": "A_HUMAN", "P2": "B_HUMAN"})
    assert table["GO term"].iloc[0] == "kinase activity"
    assert table["Proteins"].iloc[0] == "A_HUMAN, B_HUMAN"
    bars = go_terms.term_bars(result["categories"]["BP"], 5)
    assert bars["term"].iloc[0] == "kinase activity" and bars["score"].iloc[0] >= 0


def test_go_lookup_failure_is_reported(tmp_path):
    def fetch(_):
        raise ConnectionError("offline")

    index, error = go_terms.term_index(tmp_path, ["P1"], fetch=fetch)
    assert index is None and "offline" in error


def test_short_sample_names():
    from src.common.samples import short_sample_names

    names = ["B1_Ap_DIA_ratio01_DDA", "B2_Ap_DIA_ratio01_DDA_2", "B3_Ap_DIA_ratio10_DDA"]
    assert list(short_sample_names(names).values()) == ["B1_ratio01", "B2_ratio01_2", "B3_ratio10"]
    assert short_sample_names(["a_x", "a_x"]) == {"a_x": "a_x"}
    assert short_sample_names(["only_one"]) == {"only_one": "only_one"}
