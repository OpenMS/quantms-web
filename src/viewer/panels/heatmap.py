"""The most variable proteins as z-scores, in sample order (the template's
library panel), with a slider for how many."""

import streamlit as st

from src.common.insight_tables import _default_cache_dir
from src.common.guide import ensure_group_comparison
from src.common.results_helpers import (
    get_abundance_data,
    get_active_table,
    get_id_column,
    get_sample_group_map,
    log2_matrix,
)
from src.viz.heatmap import show_zscore_heatmap, top_setting
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
    head = heading("Protein Abundance Heatmap (Z-score)")
    plot = show_zscore_heatmap(
        expr_df, id_col, "quantms_protein_heatmap", cache_path=_default_cache_dir(), key=scoped("heatmap"),
        state_manager=shared_state(), groups=sample_group_map or None,
        top=top_setting(label="Proteins shown (most variable first)"),
        title=head.get("title"), help=head.get("help"), x_label="Samples", y_label="Proteins",
    )
    if plot is None:
        st.warning("Insufficient data to generate the heatmap.")
