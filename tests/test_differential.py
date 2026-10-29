"""Tests for the one-go group comparison in src/common/differential.py."""
import numpy as np
import pandas as pd
import pytest

from src.common.differential import (
    count_significant,
    default_settings,
    describe_settings,
    resolve_settings,
    run_differential,
    valid_test,
)


def _table(n_groups=2, n_proteins=60, seed=0):
    rng = np.random.default_rng(seed)
    samples = [f"sample{i + 1}[file{i + 1}]" for i in range(3 * n_groups)]
    group_map = {s: f"g{i // 3}" for i, s in enumerate(samples)}
    base = rng.lognormal(mean=20, sigma=1, size=(n_proteins, 1))
    data = base * rng.lognormal(mean=0, sigma=0.15, size=(n_proteins, len(samples)))
    data[:10, 3:6] *= 8  # ten proteins clearly up in the second group
    data[50:, :3] = 0  # ten proteins never seen in the first group
    df = pd.DataFrame(data, columns=samples)
    df.insert(0, "ProteinName", [f"P{i:05d}" for i in range(n_proteins)])
    return df, group_map


def test_defaults_run_end_to_end_and_find_the_changed_proteins():
    df, groups = _table()
    out = run_differential(df, groups, "ProteinName", default_settings())
    stats = out["statistics_df"]
    assert {"log2FC", "p-value", "p-adj"} <= set(stats.columns)
    assert out["normalized_transform"] == "log2"
    up, _ = count_significant(stats)
    assert up >= 5
    top = stats.sort_values("p-value")["ProteinName"].head(10)
    assert top.isin([f"P{i:05d}" for i in range(10)]).sum() >= 5


def test_continues_from_a_step_done_by_hand():
    df, groups = _table()
    filtered = df.head(30)
    out = run_differential(df, groups, "ProteinName", default_settings(), done={"filtered_df": filtered})
    assert out["filtered_df"] is filtered
    assert len(out["statistics_df"]) == 30


def test_needs_two_groups():
    df, groups = _table()
    with pytest.raises(ValueError, match="two groups"):
        run_differential(df, {s: "g0" for s in groups}, "ProteinName", default_settings())


def test_test_falls_back_to_limma_for_more_groups():
    assert valid_test("welch", 3) == "limma_like"
    assert valid_test("anova", 3) == "anova"
    assert valid_test("anova", 2) == "limma_like"
    df, groups = _table(n_groups=3)
    settings = {**default_settings(), "postproc-stat-test": "welch"}
    out = run_differential(df, groups, "ProteinName", settings)
    assert out["stat_test"] == "limma_like"


def test_saved_settings_win_over_defaults_and_are_described():
    settings = resolve_settings({"postproc-fdr": "Bonferroni"})
    assert settings["postproc-fdr"] == "Bonferroni"
    assert settings["postproc-transform"] == "log2"
    lines = describe_settings(settings)
    assert len(lines) == 4 and "Bonferroni" in lines[3]


def test_comparison_direction_matches_the_sign_of_log2fc():
    from src.common.differential import comparison_groups

    df, groups = _table()
    ref, comp = comparison_groups(groups, list(groups))
    assert (ref, comp) == ("g0", "g1")
    stats = run_differential(df, groups, "ProteinName", default_settings())["statistics_df"]
    # The planted proteins are higher in g1, so they must come out positive
    assert (stats.set_index("ProteinName").loc[[f"P{i:05d}" for i in range(10)], "log2FC"] > 0).all()
