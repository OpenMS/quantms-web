"""Variation within groups: a box per group of its proteins' CVs (the
library's replicate boxes), the picked protein marked on every box."""

import polars as pl
import streamlit as st
from openms_insight import SeriesPlot

from src.common.differential import lead_accession
from src.common.insight_tables import _default_cache_dir
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from src.viewer.cv import LOW_CV, MIN_VALUES, cv_summary, group_cvs, plot_table
from src.viewer.scope import heading, needs_results, scoped, shared_state, zoom
from src.viz.plumbing import group_colors
from src.viz.quantification import replicate_cv
from src.viz.series import show_series

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
        # Samples without groups: one box over all of them
        groups = {s: ALL_SAMPLES for s in sample_cols}
    mat = pivot_df.set_index(id_col)[sample_cols]
    cvs = group_cvs(mat, groups)
    if cvs.empty:
        st.info(f"The variation within a group needs a group of at least {MIN_VALUES} samples.")
        return
    summary = cv_summary(cvs, groups)
    head = heading()
    settings = replicate_cv.settings(
        pl.from_pandas(plot_table(cvs, lead_accession)), group_colors(groups.values()), selection="protein",
        id_column="accession", x_label="Group",
    )
    if head.get("title"):
        settings["title"] = settings["title"].replace("Variation between replicates (CV)", head["title"])
    data = settings.pop("data")
    show_series(
        data, "quantms_cv_distribution", cache_path=_default_cache_dir(), key=scoped("cv_distribution"),
        state_manager=shared_state(), height=400, help=head.get("help"),
        zoom=zoom(SeriesPlot, "cv_distribution").get("zoom_identifier"), **settings,
    )
    lines = [
        f"**{row.group}** ({row.samples} samples): median CV {row.median:.1f} %, "
        f"{row.below:.0f} % of {row.proteins:,} proteins below {LOW_CV:.0f} %"
        for row in summary.itertuples()
    ]
    st.markdown("  \n".join(lines))
