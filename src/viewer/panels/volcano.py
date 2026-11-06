"""results_volcano as a Viewer component."""

import streamlit as st

import polars as pl
from src.common.differential import comparison_groups, raw_p_cutoff
from src.common.guide import ensure_group_comparison, show_comparison_settings
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from src.common.insight_tables import plot_cache
from openms_insight import VolcanoPlot
from src.viewer.scope import needs_results, scoped, shared_state


def render(ctx) -> None:
    params = ctx.params
    result = get_abundance_data(ctx.workspace)
    if result is None:
        return needs_results()
    # Computed with the recommended settings unless Edit comparison already did
    statistics_df = ensure_group_comparison(params)
    if statistics_df is None:
        return
    if statistics_df.empty:
        st.info("No data available for volcano plot.")
        return
    pivot_df, expr_df, group_map = result
    id_col = get_id_column(ctx.workspace, pivot_df)
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    direction = comparison_groups(sample_group_map, [c for c in pivot_df.columns if c in sample_group_map])
    if direction:
        st.markdown(f"Right: higher in **{direction[1]}**. Left: higher in **{direction[0]}**.")
    # 2. Clean data and convert to Polars for component input
    volcano_df = statistics_df.dropna(subset=["log2FC", "p-value", "p-adj"]).copy()
    # Short labels: "sp|P05161|ISG15_HUMAN" -> "ISG15_HUMAN", so neighbouring labels overlap less
    volcano_df["label"] = volcano_df[id_col].astype(str).str.split("|").str[-1]
    volcano_pl_lazy = pl.from_pandas(volcano_df).lazy()
    # 3. Configure UI sliders (changing thresholds does not invalidate cache)
    fc_thresh = st.slider(
        "Minimum change (|log2 fold change|)",
        min_value=0.5,
        max_value=3.0,
        value=1.0,
        step=0.1,
        help="1 = at least twice (or half) as abundant.",
        key=scoped("volcano_fc"),
    )
    p_thresh = st.slider(
        "Maximum adjusted p-value",
        min_value=0.001,
        max_value=0.1,
        value=0.05,
        step=0.001,
        help="0.05 means about 5 % of the proteins called significant are expected to be false positives.",
        key=scoped("volcano_p"),
    )
    # 4. Initialize the OpenMS-Insight VolcanoPlot component
    volcano_plot_component = VolcanoPlot(
        **plot_cache("quantms_volcano_plot", volcano_df),
        data=volcano_pl_lazy,
        log2fc_column="log2FC",
        pvalue_column="p-value",
        label_column="label",
        y_label="-log10(p-value)",
        up_color="#E74C3C",
        down_color="#3498DB",
        ns_color="#95A5A6",
        show_threshold_lines=True,
        threshold_line_style="dash",
    )
    # 5. Render the component
    state_manager = shared_state()
    volcano_plot_component(
        state_manager=state_manager,
        fc_threshold=fc_thresh,
        # The adjusted-p cut-off as a raw p-value, since the y-axis is the raw p-value
        p_threshold=raw_p_cutoff(volcano_df, p_thresh),
        max_labels=10,  # Display labels for the top N significant proteins
        height=600,
        key=scoped("volcano_plot"),
    )
    up_count = ((volcano_df["p-adj"] <= p_thresh) & (volcano_df["log2FC"] >= fc_thresh)).sum()
    down_count = ((volcano_df["p-adj"] <= p_thresh) & (volcano_df["log2FC"] <= -fc_thresh)).sum()
    if direction:
        st.markdown(f"**Higher in {direction[1]}:** {up_count} | **Lower in {direction[1]}:** {down_count}")
    else:
        st.markdown(f"**Up-regulated:** {up_count} | **Down-regulated:** {down_count}")
    show_comparison_settings()
