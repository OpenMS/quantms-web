"""Abundance of the picked protein in every sample, coloured by group."""
import numpy as np
import pandas as pd
import polars as pl
import streamlit as st
from openms_insight import SeriesPlot

from src.common.differential import lead_accession
from src.common.insight_tables import plot_cache
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from src.viewer.scope import needs_results, scoped, shared_state, zoom


def profile_table(pivot_df: pd.DataFrame, id_col: str, sample_group_map: dict) -> pd.DataFrame:
    """One row per protein and sample: sample number, log2 abundance, group, for a SeriesPlot."""
    samples = [c for c in pivot_df.columns if c in sample_group_map]
    long = pivot_df[[id_col, *samples]].melt(id_vars=id_col, var_name="sample", value_name="abundance")
    long["abundance"] = pd.to_numeric(long["abundance"], errors="coerce")
    long = long[long["abundance"] > 0]
    order = {s: i + 1 for i, s in enumerate(samples)}
    return pd.DataFrame({
        "accession": long[id_col].map(lead_accession),
        "x": long["sample"].map(order).astype(float),
        "y": np.log2(long["abundance"].astype(float)),
        "group": long["sample"].map(sample_group_map).astype(str),
        "sample": long["sample"].astype(str),
    }).sort_values(["accession", "x"]).reset_index(drop=True)


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
    if not state_manager.get_selection("protein"):
        st.caption("Pick a protein in the volcano plot or a table to see its abundance in each sample.")
    samples = table.drop_duplicates("x").sort_values("x")
    st.caption("Samples: " + " · ".join(f"{int(x)} {s}" for x, s in zip(samples["x"], samples["sample"])))
    SeriesPlot(
        **plot_cache("quantms_abundance_profile", table),
        **zoom(SeriesPlot, "abundance_profile"),
        data=pl.from_pandas(table).lazy(),
        filters={"protein": "accession"},
        # Nothing until a protein is picked, rather than every protein at once
        filter_defaults={"protein": "__none__"},
        x_column="x",
        y_column="y",
        series_column="group",
        hover_column="sample",
        mode="markers",
        x_label="Sample",
        y_label="log2 abundance",
        y_range_mode="normal",
    )(key=scoped("abundance_profile"), state_manager=state_manager, height=360)
