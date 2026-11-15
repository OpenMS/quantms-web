"""results_pca as a Viewer component."""

import streamlit as st

from src.common.guide import ensure_group_comparison
from src.common.insight_tables import _default_cache_dir
from src.common.results_helpers import (
    get_abundance_data,
    get_active_table,
    get_id_column,
    get_sample_group_map,
    sample_columns,
)
from openms_insight import PCAPlot
from openms_insight.core.settings import setting_value
from src.viz.pca import pca_component as build_pca
from src.viewer.scope import heading, needs_results, scoped, shared_state, zoom


def render(ctx) -> None:
    params = ctx.params
    # 1. Load abundance data (base wide-format table + sample -> group mapping)
    result = get_abundance_data(ctx.workspace)
    if result is None:
        return needs_results()
    pivot_df, expr_df, group_map = result
    id_col = get_id_column(ctx.workspace, pivot_df)
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    # --- STEP 1: Use the latest downstream step's output (same table as Statistics) ---
    # Plots follow the group comparison, run with the recommended settings if needed
    ensure_group_comparison(params, show_errors=False)
    base_df, _, _ = get_active_table(pivot_df)
    # 2. Extract active sample columns and detect unique biological groups
    sample_cols = sample_columns(base_df, sample_group_map)
    unique_groups = sorted({sample_group_map[s] for s in sample_cols if s in sample_group_map})
    if len(sample_cols) < 2:
        st.info("PCA requires at least 2 samples.")
        return
    if len(unique_groups) < 2:
        st.warning(
            "Only one biological group was detected - points will still be plotted, "
            "but group-based coloring requires 2 or more groups."
        )
    expr_df_wide = base_df.set_index(id_col)[sample_cols]
    max_available = expr_df_wide.shape[0]
    state_manager = shared_state()
    # The plot's settings live in its toolbar (⚙): proteins used and the axes
    settings = []
    top_n = max_available
    if max_available > 20:
        top_setting = {
            "identifier": "pca_top_n", "label": "Proteins used (most variable first)", "kind": "range",
            "min": 20, "max": min(5000, max_available), "step": 10, "default": min(500, max_available),
        }
        settings.append(top_setting)
        top_n = int(setting_value(state_manager, top_setting))
    top_proteins = expr_df_wide.var(axis=1).sort_values(ascending=False).head(top_n).index
    expr_df_pca = expr_df_wide.loc[top_proteins].reset_index()
    if expr_df_pca.shape[0] < 2:
        st.info("Not enough proteins after variance filtering for PCA.")
        return
    # 3. The library's PCA panel (computes PCA internally, samples coloured by group)
    head = heading("Sample PCA")
    try:
        pca_component = build_pca(
            expr_df_pca, sample_group_map, "quantms_pca_plot", cache_path=_default_cache_dir(),
            samples=sample_cols, zoom=zoom(PCAPlot, "pca").get("zoom_identifier"), n_components=5,
            title=head.get("title") or "Sample PCA",
        )
    except ValueError as e:
        st.error(f"PCA computation failed: {e}")
        return
    pc_columns = pca_component.get_pc_columns()
    # Which component pair to view (no recomputation needed)
    x_setting = {"identifier": "pca_pc_x", "label": "X axis", "kind": "choice", "options": pc_columns,
                 "default": pc_columns[0]}
    y_setting = {"identifier": "pca_pc_y", "label": "Y axis", "kind": "choice", "options": pc_columns,
                 "default": pc_columns[1 if len(pc_columns) > 1 else 0]}
    settings += [x_setting, y_setting]
    pc_x = int(str(setting_value(state_manager, x_setting)).replace("PC", ""))
    pc_y = int(str(setting_value(state_manager, y_setting)).replace("PC", ""))
    pca_component.set_toolbar_settings(settings)
    pca_component(state_manager=state_manager, pc_x=pc_x, pc_y=pc_y, height=520, key=scoped("pca_plot"),
                  help=head.get("help"))

