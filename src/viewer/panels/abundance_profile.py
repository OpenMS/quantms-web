"""Abundance of the picked protein in every sample, coloured by group."""
import inspect

import numpy as np
import pandas as pd
import polars as pl
import streamlit as st
from openms_insight import SeriesPlot

from src.common.differential import lead_accession
from src.common.insight_tables import plot_cache
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from src.viewer.scope import heading, needs_results, scoped, shared_state, zoom


# The colours the volcano and PCA use for the two groups first
PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]


def profile_table(pivot_df: pd.DataFrame, id_col: str, sample_group_map: dict) -> pd.DataFrame:
    """One row per protein and sample (sample number, log2 abundance, group), for a SeriesPlot.

    The samples are numbered group by group, so each group's samples sit side by
    side, and each group gets a "<group> mean" line across its samples: the gap
    between the two means is the fold change the volcano plot shows.
    """
    samples = [c for c in pivot_df.columns if c in sample_group_map]
    groups = list(dict.fromkeys(str(sample_group_map[s]) for s in samples))
    samples = sorted(samples, key=lambda s: groups.index(str(sample_group_map[s])))
    long = pivot_df[[id_col, *samples]].melt(id_vars=id_col, var_name="sample", value_name="abundance")
    long["abundance"] = pd.to_numeric(long["abundance"], errors="coerce")
    long = long[long["abundance"] > 0]
    order = {s: i + 1 for i, s in enumerate(samples)}
    points = pd.DataFrame({
        "accession": long[id_col].map(lead_accession),
        "x": long["sample"].map(order).astype(float),
        "y": np.log2(long["abundance"].astype(float)),
        "group": long["sample"].map(sample_group_map).astype(str),
        "sample": long["sample"].astype(str),
    })
    if points.empty:
        return points.reset_index(drop=True)
    span = {g: (min(order[s] for s in samples if str(sample_group_map[s]) == g),
                max(order[s] for s in samples if str(sample_group_map[s]) == g)) for g in groups}
    means = points.groupby(["accession", "group"], as_index=False)["y"].mean()
    lines = pd.concat([
        means.assign(x=means["group"].map(lambda g, i=i: float(span[g][i])))
        for i in (0, 1)
    ])
    lines["sample"] = "mean of " + lines["group"]
    lines["group"] = lines["group"] + " mean"
    return (pd.concat([points, lines[points.columns]])
            .sort_values(["accession", "group", "x"]).reset_index(drop=True))


def render(ctx) -> None:
    result = get_abundance_data(ctx.workspace)
    if result is None:
        return needs_results()
    pivot_df, _, group_map = result
    id_col = get_id_column(ctx.workspace, pivot_df)
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    table = profile_table(pivot_df, id_col, sample_group_map)
    if table.empty:
        st.info("No abundances in this run.")
        return
    state_manager = shared_state()
    # A group and its mean line share one colour
    groups = [g for g in dict.fromkeys(table["group"]) if not g.endswith(" mean")]
    colors = {}
    for i, group in enumerate(groups):
        colors[group] = colors[f"{group} mean"] = PALETTE[i % len(PALETTE)]
    # Samples are points, not a sequence: only the group means are lines
    modes = {}
    if "series_modes" in inspect.signature(SeriesPlot.__init__).parameters:
        modes = {"series_modes": {f"{group} mean": "lines" for group in groups}}
    SeriesPlot(
        **plot_cache("quantms_abundance_profile_v2", table),
        **zoom(SeriesPlot, "abundance_profile"),
        **heading(),
        data=pl.from_pandas(table).lazy(),
        filters={"protein": "accession"},
        # Nothing until a protein is picked, rather than every protein at once
        filter_defaults={"protein": "__none__"},
        x_column="x",
        y_column="y",
        series_column="group",
        series_colors=colors,
        hover_column="sample",
        mode="markers",
        **modes,
        x_label="Sample",
        y_label="log2 abundance",
        y_range_mode="normal",
    )(key=scoped("abundance_profile"), state_manager=state_manager, height=560)
