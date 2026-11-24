"""results_heatmap_clustered as a Viewer component."""

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
from openms_insight.core.settings import setting_value
from src.viz.heatmap import clustered_height, show_clustered_heatmap, zscores
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
    heatmap_z = zscores(expr_df, top=top_n)
    if heatmap_z.empty:
        st.warning("Insufficient data to generate the heatmap.")
        return
    head = heading("Protein Abundance Heatmap (Z-score, clustered)")
    # Height grows with the number of proteins so row labels stay readable
    show_clustered_heatmap(
        heatmap_z.reset_index(), sample_group_map, id_col, "quantms_clustered_heatmap",
        cache_path=_default_cache_dir(), key=scoped("clustered_heatmap"), state_manager=state_manager,
        height=clustered_height(top_n), help=head.get("help"), toolbar_settings=[top_setting],
        samples=expr_df.columns.tolist(), title=head.get("title"), y_label="Proteins",
    )
