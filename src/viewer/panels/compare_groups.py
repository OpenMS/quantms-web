"""The group comparison table (the template's library panel)."""

import streamlit as st
from src.common.differential import comparison_groups, lead_accession
from src.common.guide import ensure_group_comparison, page_link
from src.common.insight_tables import _default_cache_dir
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from src.viewer.scope import GO_SET, change_cutoffs, heading, needs_results, scoped, shared_state
from src.viz.quantification.comparison_table import comparison_title, show_comparison_table


def render(ctx) -> None:
    params = ctx.params
    result = get_abundance_data(ctx.workspace)
    if result is None:
        return needs_results()
    pivot_df, _, group_map = result
    id_col = get_id_column(ctx.workspace, pivot_df)
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    sample_cols = [c for c in pivot_df.columns if c in sample_group_map]
    if any(sample_group_map[s] == "Unassigned" for s in sample_cols):
        st.warning("Samples without a group are left out. Set their group on **Describe experiment**.")
        page_link("configure")
    statistics_df = ensure_group_comparison(params)
    if statistics_df is None:
        return
    direction = comparison_groups(sample_group_map, sample_cols)
    fc, padj = change_cutoffs(params)
    reference, compared = direction or (None, None)
    table = statistics_df[[id_col, "log2FC", "p-adj", "p-value"]].copy()
    # Picking a row picks its lead protein in every component of the run
    # (the Identification tables, the volcano plot), and the other way round
    table["accession"] = table[id_col].map(lead_accession)
    head = heading(comparison_title(reference, compared), detail=comparison_title(reference, compared) if direction else None)
    show_comparison_table(
        table, id_col, "Protein", "quantms_compare_groups", cache_path=_default_cache_dir(),
        key=scoped("compare_groups"), state_manager=shared_state(), reference=reference, compared=compared,
        fc=fc, padj=padj, title=head.get("title"), help=head.get("help"),
        interactivity={"protein": "accession"},
        # The check boxes pick proteins for the GO enrichment, like the volcano plot
        multi_select={GO_SET: "accession"},
    )
