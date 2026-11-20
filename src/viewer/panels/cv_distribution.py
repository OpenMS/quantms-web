"""Variation within groups: the distribution of the proteins' CV in each group."""

import polars as pl
import streamlit as st

from openms_insight import SeriesPlot
from src.common.cv import LOW_CV, MIN_VALUES, cv_distribution, cv_summary, group_cvs
from src.common.insight_tables import plot_cache
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from src.viewer.panels.heatmap_clustered import group_colors
from src.viewer.scope import needs_results, scoped, shared_state, zoom

ALL_SAMPLES = "All samples"


def render(ctx) -> None:
    result = get_abundance_data(ctx.workspace)
    if result is None:
        return needs_results()
    pivot_df, expr_df, group_map = result
    id_col = get_id_column(ctx.workspace, pivot_df)
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    # The measured abundances, not the group comparison's table: filled-in
    # values would make the replicates look more alike than they are.
    sample_cols = [c for c in expr_df.columns if c in pivot_df.columns]
    groups = {s: sample_group_map[s] for s in sample_cols if s in sample_group_map}
    if not groups:
        # Samples without groups: one distribution over all of them
        groups = {s: ALL_SAMPLES for s in sample_cols}
    mat = pivot_df.set_index(id_col)[sample_cols]
    cvs = group_cvs(mat, groups)
    if cvs.empty:
        st.info(f"The variation within a group needs a group of at least {MIN_VALUES} samples.")
        return
    summary = cv_summary(cvs, groups)
    st.caption(
        "Coefficient of variation (standard deviation divided by the mean, in %) of each protein's "
        f"measured abundances in the samples of a group, over the proteins measured in at least {MIN_VALUES} "
        "of them. Replicates measured precisely give most proteins a small CV."
    )
    curves = cv_distribution(cvs)
    colors = group_colors(groups)
    plot = SeriesPlot(
        **plot_cache("quantms_cv_distribution", curves),
        **zoom(SeriesPlot, "cv_distribution"),
        data=pl.from_pandas(curves).lazy(),
        x_column="cv",
        y_column="proteins",
        series_column="group",
        series_colors={g: colors[g] for g in summary["group"]},
        mode="lines+markers",
        title="CV of the proteins in each group",
        x_label="CV (%)",
        y_label="Proteins (%)",
    )
    plot(key=scoped("cv_distribution"), state_manager=shared_state(), height=380)
    lines = [
        f"**{row.group}** ({row.samples} samples): median CV {row.median:.1f} %, "
        f"{row.below:.0f} % of {row.proteins:,} proteins below {LOW_CV:.0f} %"
        for row in summary.itertuples()
    ]
    st.markdown("  \n".join(lines))
