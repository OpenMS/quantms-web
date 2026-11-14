"""results_pca as a Viewer component."""

import polars as pl
import streamlit as st

from src.common.guide import ensure_group_comparison
from src.common.insight_tables import plot_cache
from src.common.results_helpers import (
    get_abundance_data,
    get_active_table,
    get_id_column,
    get_sample_group_map,
    sample_columns,
)
from openms_insight import PCAPlot
from openms_insight.core.settings import setting_value
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
    # Prepare structural Polars metadata DataFrame required by PCAPlot
    metadata_pl = pl.DataFrame(
        [{"sample_id": s, "group": sample_group_map[s]} for s in sample_cols if s in sample_group_map],
        schema={"sample_id": pl.String, "group": pl.String},
    )
    pca_lazy = pl.from_pandas(expr_df_pca).lazy()
    # 3. Initialize the OpenMS-Insight PCAPlot component (computes PCA internally)
    try:
        pca_component = PCAPlot(
            **plot_cache("quantms_pca_plot", expr_df_pca),
            **zoom(PCAPlot, "pca"),
            data=pca_lazy,
            metadata=metadata_pl,
            sample_id_field="sample_id",
            group_field="group",
            n_components=5,
            **heading("Sample PCA"),
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
    pca_component(state_manager=state_manager, pc_x=pc_x, pc_y=pc_y, height=520, key=scoped("pca_plot"))

