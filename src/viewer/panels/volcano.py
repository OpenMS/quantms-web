"""results_volcano as a Viewer component."""

import streamlit as st
from openms_insight import VolcanoPlot

from src.common.differential import comparison_groups, lead_accession, raw_p_cutoff
from src.common.guide import ensure_group_comparison
from src.common.insight_tables import _default_cache_dir
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from src.viz.quantification.volcano import short_label, show_volcano
from src.viewer.scope import GO_SET, heading, change_cutoffs, needs_results, scoped, shared_state, zoom


def render(ctx) -> None:
    result = get_abundance_data(ctx.workspace)
    if result is None:
        return needs_results()
    # Computed with the recommended settings unless Edit comparison already did
    statistics_df = ensure_group_comparison(ctx.params)
    if statistics_df is None:
        return
    if statistics_df.empty:
        st.info("No data available for volcano plot.")
        return
    pivot_df, _, group_map = result
    id_col = get_id_column(ctx.workspace, pivot_df)
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    direction = comparison_groups(sample_group_map, [c for c in pivot_df.columns if c in sample_group_map])
    fc, padj = change_cutoffs(ctx.params)
    volcano_df = statistics_df.dropna(subset=["log2FC", "p-value", "p-adj"]).copy()
    # Short labels: "sp|P05161|ISG15_HUMAN" -> "ISG15_HUMAN", so neighbouring labels overlap less
    volcano_df["label"] = volcano_df[id_col].map(short_label)
    # A clicked point picks its lead protein in every component of the run
    volcano_df["accession"] = volcano_df[id_col].map(lead_accession)
    reference, compared = direction or (None, None)
    head = heading(f"{compared} vs. {reference}" if direction else None,
                   detail=f"{compared} vs. {reference}" if direction else None)
    show_volcano(
        volcano_df, "quantms_volcano_v2", cache_path=_default_cache_dir(), key=scoped("volcano_plot"),
        state_manager=shared_state(), zoom=zoom(VolcanoPlot, "volcano").get("zoom_identifier"),
        fc_threshold=fc,
        # The adjusted-p cut-off as a raw p-value, since the y-axis is the raw p-value
        p_threshold=raw_p_cutoff(volcano_df, padj),
        help=head.get("help"), title=head.get("title"), reference=reference, compared=compared,
        interactivity={"protein": "accession"},
        # All changed / Higher / Lower above the plot, or a box or lasso, picks
        # proteins for the GO enrichment
        multi_select={GO_SET: "accession"},
    )
