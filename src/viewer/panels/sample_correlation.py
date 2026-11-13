"""Sample correlation: how alike the samples' protein abundances are, as a clustered heatmap."""

import polars as pl
import streamlit as st

from openms_insight import ClusteredHeatmap
from src.common.insight_tables import plot_cache
from src.common.results_helpers import (
    get_abundance_data,
    get_id_column,
    get_sample_group_map,
    log2_matrix,
)
from src.common.sample_correlation import (
    least_similar,
    median_correlation,
    median_shared_proteins,
    sample_correlation,
    shared_proteins,
)
from src.viewer.panels.heatmap_clustered import group_colors
from src.viewer.scope import needs_results, scoped, shared_state

SAMPLE_COL = "sample"


def render(ctx) -> None:
    result = get_abundance_data(ctx.workspace)
    if result is None:
        return needs_results()
    pivot_df, expr_df, group_map = result
    id_col = get_id_column(ctx.workspace, pivot_df)
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    # The measured abundances, not the group comparison's table: values filled
    # in for missing ones would make samples look more alike than they are.
    sample_cols = [c for c in expr_df.columns if c in pivot_df.columns]
    if len(sample_cols) < 2:
        st.info("Sample correlation needs at least two samples.")
        return
    mat = log2_matrix(pivot_df, id_col, sample_cols, is_log2=False)
    corr = sample_correlation(mat)
    if corr.isna().all().all():
        st.info("No two samples share enough measured proteins to compare them.")
        return
    shared = shared_proteins(mat)
    st.caption(
        "Pearson correlation of the log2 protein abundances of every two samples, over the proteins "
        f"measured in both (median {median_shared_proteins(shared):,} proteins per pair). "
        "Replicates of one group should be the most alike; samples are ordered by similarity."
    )
    corr_df = corr.rename_axis(SAMPLE_COL).reset_index()
    metadata = pl.DataFrame(
        [{"sample_id": s, "group": sample_group_map[s]} for s in corr.columns if s in sample_group_map],
        schema={"sample_id": pl.String, "group": pl.String},
    )
    component = ClusteredHeatmap(
        **plot_cache("quantms_sample_correlation", corr_df),
        id_col=SAMPLE_COL,
        data=pl.from_pandas(corr_df).lazy(),
        metadata=metadata,
        row_cluster=True,
        col_cluster=True,
        title="Sample correlation (Pearson, log2 abundances)",
        x_label="Samples",
        y_label="Samples",
        intensity_label="Pearson r",
        colorscale="Viridis",
        reversescale=False,
        group_colors=group_colors(sample_group_map),
    )
    height = max(500, min(1200, 250 + 28 * len(corr)))
    component(state_manager=shared_state(), height=height, key=scoped("sample_correlation"))
    medians = median_correlation(corr).dropna()
    outliers = least_similar(corr)
    if outliers:
        st.warning(
            "Less alike than the other samples: "
            + ", ".join(f"**{s}** (median r {medians[s]:.2f})" for s in outliers)
            + f", against a median of {medians.median():.2f}. Check whether these runs worked and belong to their group."
        )
    elif len(medians):
        st.caption(f"Median correlation between samples: {medians.median():.2f}, lowest {medians.iloc[0]:.2f} ({medians.index[0]}).")

