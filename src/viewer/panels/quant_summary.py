"""Quantification summary: proteins, samples, groups, proteins up and down, missing values as tiles."""
import pandas as pd
import streamlit as st

from src.common.differential import comparison_groups, count_significant
from src.common.guide import ensure_group_comparison
from src.common.results_helpers import get_abundance_data, get_sample_group_map
from src.viewer.scope import change_cutoffs, needs_results


def render(ctx) -> None:
    result = get_abundance_data(ctx.workspace)
    if result is None:
        return needs_results()
    pivot_df, expr_df, group_map = result
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    sample_cols = [c for c in pivot_df.columns if c in sample_group_map]
    groups = {g for g in sample_group_map.values() if str(g).strip() and g != "Unassigned"}
    # Missing in the report (expr_df keeps only proteins measured in every sample)
    measured = [c for c in (sample_cols or list(expr_df.columns)) if c in pivot_df.columns]
    values = pivot_df[measured].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float) if measured else None
    missing = None
    if values is not None and values.size:
        missing = round(100 * float(((values != values) | (values == 0)).mean()), 1)
    up = down = None
    direction = comparison_groups(sample_group_map, sample_cols) if len(groups) >= 2 else None
    fc, padj = change_cutoffs(ctx.params)
    if direction:
        statistics_df = ensure_group_comparison(ctx.params, show_errors=False)
        if statistics_df is not None and not statistics_df.empty:
            up, down = count_significant(statistics_df, fc, padj)
        st.caption(
            f"Up and down: higher or lower in {direction[1]} than in {direction[0]}, "
            f"|log2 fold change| ≥ {fc:g} and adjusted p-value ≤ {padj:g} (Edit comparison, step 4)."
        )

    def show(value, unit: str = "") -> str:
        return "–" if value is None else f"{value:,}{unit}"

    tiles = [
        ("Quantified proteins", show(len(pivot_df)), "Without contaminants."),
        ("Samples", show(len(measured) or len(expr_df.columns)), None),
        ("Groups", show(len(groups) or None), "Set on Describe experiment."),
        ("Proteins up", show(up), None),
        ("Proteins down", show(down), None),
        ("Missing values", show(missing, " %"), "Share of protein abundances not measured, before imputation."),
    ]
    for col, (label, value, help_text) in zip(st.columns(len(tiles)), tiles):
        col.metric(label, value, help=help_text)
