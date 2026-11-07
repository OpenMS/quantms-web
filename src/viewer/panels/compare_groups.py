"""compare_groups as a Viewer component."""

import streamlit as st
from src.common.differential import comparison_groups, count_significant, lead_accession
from src.common.guide import (
    ensure_group_comparison,
    page_link,
    show_comparison_settings,
)
from src.common.insight_tables import show_insight_table
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from src.viewer.scope import needs_results, shared_state


def render(ctx) -> None:
    params = ctx.params
    result = get_abundance_data(ctx.workspace)
    if result is None:
        return needs_results()
    pivot_df, _, group_map = result
    id_col = get_id_column(ctx.workspace, pivot_df)
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    sample_cols = [c for c in pivot_df.columns if c in sample_group_map]
    st.subheader("Your groups")
    groups: dict[str, list[str]] = {}
    for s in sample_cols:
        groups.setdefault(sample_group_map[s], []).append(s)
    cols = st.columns(max(len(groups), 1))
    for col, (group, members) in zip(cols, sorted(groups.items())):
        col.metric(group, f"{len(members)} sample{'' if len(members) == 1 else 's'}")
    if "Unassigned" in groups:
        st.warning("Samples without a group are left out. Set their group on **Describe experiment**.")
        page_link("configure")
    statistics_df = ensure_group_comparison(params)
    if statistics_df is None:
        return
    direction = comparison_groups(sample_group_map, sample_cols)
    up, down = count_significant(statistics_df)
    st.subheader("Result")
    m = st.columns(3)
    m[0].metric("Proteins compared", f"{statistics_df['p-adj'].notna().sum():,}")
    if direction:
        ref, comp = direction
        m[1].metric(f"Higher in {comp}", f"{up:,}", help="Adjusted p-value ≤ 0.05 and at least twice as abundant.")
        m[2].metric(f"Lower in {comp}", f"{down:,}", help="Adjusted p-value ≤ 0.05 and at most half as abundant.")
        note = f"A positive log2 fold change means higher in **{comp}** than in **{ref}**."
        if len(groups) > 2:
            note += (
                " With more than two groups the p-value tests whether any group differs; the fold "
                f"change still compares **{comp}** with **{ref}**."
            )
        st.caption(note + " Counts: |log2 fold change| ≥ 1, adjusted p-value ≤ 0.05.")
    show_comparison_settings()
    st.subheader("All compared proteins, most significant first")
    table = statistics_df.sort_values("p-adj", na_position="last")[[id_col, "log2FC", "p-adj", "p-value"]].copy()
    if direction:
        sig = (table["p-adj"] <= 0.05) & (table["log2FC"].abs() >= 1)
        table.insert(1, "Change", "–")
        table.loc[sig & (table["log2FC"] > 0), "Change"] = f"higher in {direction[1]}"
        table.loc[sig & (table["log2FC"] < 0), "Change"] = f"lower in {direction[1]}"
    for col in ("p-adj", "p-value"):
        # As text with three significant digits: the table stores floats as float32,
        # which would show 0.155 as 0.1550000011920929. Rows are already in p-value order.
        table[col] = table[col].map(lambda v: "" if v != v else f"{v:.3g}")
    table = table.rename(
        columns={
            id_col: "Protein",
            "log2FC": "log2 fold change",
            "p-adj": "adjusted p-value",
        }
    )
    # Picking a row picks its lead protein in every component of the run
    # (the Identification tables, the volcano plot), and the other way round
    table["accession"] = table["Protein"].map(lead_accession)
    show_insight_table(
        table,
        name="compare-groups",
        layout="fitColumns",
        state_manager=shared_state(),
        interactivity={"protein": "accession"},
        column_overrides={
            "log2 fold change": {"formatter": "money", "formatterParams": {"thousand": False, "precision": 2}},
            "accession": {"visible": False},
        },
    )
