import numpy as np
import pandas as pd

from src.common.sample_correlation import (
    least_similar,
    median_correlation,
    median_shared_proteins,
    sample_correlation,
    shared_proteins,
)


def _matrix(seed=0, n=200):
    rng = np.random.default_rng(seed)
    base = rng.normal(25, 2, n)
    cols = {}
    for name in ["ctrl_1", "ctrl_2", "ctrl_3", "drug_1", "drug_2", "drug_3"]:
        cols[name] = base + rng.normal(0, 0.2, n)
    return pd.DataFrame(cols, index=[f"P{i}" for i in range(n)])


def test_replicates_correlate():
    corr = sample_correlation(_matrix())
    assert list(corr.columns) == list(corr.index)
    assert np.allclose(np.diag(corr), 1.0)
    assert corr.min().min() > 0.95
    assert least_similar(corr) == []


def test_pairwise_complete_and_min_shared():
    mat = _matrix()
    mat.loc[mat.index[:195], "drug_3"] = np.nan  # shares only 5 proteins with the others
    corr = sample_correlation(mat)
    assert corr["drug_3"].drop("drug_3").isna().all()
    assert corr.loc["ctrl_1", "ctrl_2"] > 0.95
    shared = shared_proteins(mat)
    assert shared.loc["ctrl_1", "drug_3"] == 5
    assert shared.loc["ctrl_1", "ctrl_2"] == 200


def test_infinite_values_are_missing():
    mat = _matrix()
    mat.iloc[0, 0] = -np.inf  # log2 of a zero intensity
    assert np.isfinite(sample_correlation(mat).loc["ctrl_1", "ctrl_2"])


def test_odd_sample_is_named():
    mat = _matrix()
    rng = np.random.default_rng(1)
    mat["drug_3"] = rng.normal(25, 2, len(mat))  # unrelated to the others
    corr = sample_correlation(mat)
    assert least_similar(corr) == ["drug_3"]
    assert median_correlation(corr).index[0] == "drug_3"


def test_median_shared_proteins():
    shared = pd.DataFrame([[10, 4, 6], [4, 9, 8], [6, 8, 12]], index=list("abc"), columns=list("abc"))
    assert median_shared_proteins(shared) == 6


def test_catalog_lists_component():
    from src.viewer.catalog import EXPERIMENT_WIDE, LABELS, MODULES, NEEDS

    assert LABELS["sample_correlation"] == "Sample correlation"
    assert MODULES["sample_correlation"] == "sample_correlation"
    assert NEEDS["sample_correlation"] == ["quantification", "samples"]
    assert "sample_correlation" in EXPERIMENT_WIDE
