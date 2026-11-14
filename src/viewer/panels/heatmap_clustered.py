"""results_heatmap_clustered as a Viewer component."""

import streamlit as st

import numpy as np
import polars as pl
from src.common.insight_tables import plot_cache
from src.common.guide import ensure_group_comparison
from src.common.results_helpers import (
    get_abundance_data,
    get_active_table,
    get_id_column,
    get_sample_group_map,
    log2_matrix,
)
from openms_insight import ClusteredHeatmap
from openms_insight.core.settings import setting_value
from src.viewer.scope import heading, needs_results, scoped, shared_state


def render(ctx) -> None:
    params = ctx.params
    result = get_abundance_data(ctx.workspace)
    if result is None:
        return needs_results()
    pivot_df, expr_df, group_map = result
    id_col = get_id_column(ctx.workspace, pivot_df)
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    # Plot the same protein table as PCA and Statistics: the latest downstream
    # step's output, on a log2 scale, proteins with any missing value dropped.
    # Plots follow the group comparison, run with the recommended settings if needed
    ensure_group_comparison(params, show_errors=False)
    base_df, _, is_log2 = get_active_table(pivot_df)
    # expr_df from the workflow lists exactly the sample columns in both LFQ and TMT mode
    sample_cols = [c for c in expr_df.columns if c in base_df.columns]
    expr_df = log2_matrix(base_df, id_col, sample_cols, is_log2).dropna()
    if expr_df.empty:
        st.info("No data available for heatmap.")
        return
    state_manager = shared_state()
    top_setting = {
        "identifier": "clustered_heatmap_top_n", "label": "Proteins shown (most variable first)", "kind": "range",
        "min": 10, "max": 200, "step": 10, "default": 30,
    }
    top_n = int(setting_value(state_manager, top_setting))
    var_series = expr_df.var(axis=1)
    top_proteins = var_series.sort_values(ascending=False).head(top_n).index
    heatmap_df = expr_df.loc[top_proteins]
    heatmap_z = heatmap_df.sub(heatmap_df.mean(axis=1), axis=0).div(heatmap_df.std(axis=1), axis=0)
    heatmap_z = heatmap_z.replace([np.inf, -np.inf], np.nan).dropna()
    if heatmap_z.empty:
        st.warning("Insufficient data to generate the heatmap.")
        return
    heatmap_z_reset = heatmap_z.reset_index()
    heatmap_lazy = pl.from_pandas(heatmap_z_reset).lazy()
    sample_cols = expr_df.columns.tolist()
    metadata_pl = pl.DataFrame(
        [{"sample_id": s, "group": sample_group_map[s]} for s in sample_cols if s in sample_group_map],
        schema={"sample_id": pl.String, "group": pl.String},
    )
    # Assign group annotation-bar colors in sorted-group order (matching how
    # ClusteredHeatmap._preprocess() orders unique groups internally).
    group_palette = [
        "#00BFC4",  # teal
        "#F8766D",  # salmon
        "#7CAE00",  # yellow-green
        "#C77CFF",  # lavender purple
        "#E7B800",  # gold/amber
        "#619CFF",  # blue
        "#FF61C3",  # pink/magenta
        "#00BA38",  # green
        "#FF8C42",  # orange
        "#00B0F6",  # sky blue
    ]
    unique_groups = sorted(set(sample_group_map.values()))
    group_colors = {g: group_palette[i % len(group_palette)] for i, g in enumerate(unique_groups)}
    heatmap_component = ClusteredHeatmap(
        **plot_cache("quantms_clustered_heatmap", heatmap_z),
        id_col=id_col,
        data=heatmap_lazy,
        metadata=metadata_pl,
        row_cluster=True,
        col_cluster=True,
        **heading("Protein Abundance Heatmap (Z-score, clustered)"),
        toolbar_settings=[top_setting],
        x_label="Samples",
        y_label="Proteins",
        intensity_label="Z-score",
        group_colors=group_colors,
    )
    # Scale height with the number of proteins so row labels stay readable -
    # BaseComponent otherwise defaults to a flat 400px, too short for a
    # dendrogram+heatmap composite with more than a handful of rows.
    heatmap_height = max(600, min(1400, 300 + top_n * 20))
    heatmap_component(state_manager=state_manager, height=heatmap_height, key=scoped("clustered_heatmap"))
