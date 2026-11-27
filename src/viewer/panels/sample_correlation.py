"""Sample correlation: how alike the samples' protein abundances are, as the
library's clustered sample by sample heatmap."""

import polars as pl
import streamlit as st
from openms_insight import ClusteredHeatmap

from src.common.insight_tables import _default_cache_dir
from src.common.results_helpers import (
    get_abundance_data,
    get_id_column,
    get_sample_group_map,
    log2_matrix,
)
from src.viewer.sample_correlation import (
    least_similar,
    median_correlation,
    median_shared_proteins,
    sample_correlation,
    shared_proteins,
)
from src.viewer.scope import heading, needs_results, scoped, shared_state
from src.viz.plumbing import cache_name, component, group_colors, only_accepted
from src.viz.quantification import sample_correlation as library
from src.viz.quantification.pca import sample_metadata


def correlation_source(corr) -> pl.DataFrame:
    """The library's contract: the sample in ``Sample``, one column per sample."""
    return pl.from_pandas(corr.rename_axis("Sample").reset_index())


def render(ctx) -> None:
    result = get_abundance_data(ctx.workspace)
    if result is None:
        return needs_results()
    pivot_df, expr_df, group_map = result
    id_col = get_id_column(ctx.workspace, pivot_df)
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    # The measured abundances, not the group comparison's table: values filled
    # in for missing ones would make samples look more alike than they are.
    sample_cols = [c for c in expr_df.columns if c in pivot_df.columns]
    if len(sample_cols) < 2:
        st.info("Sample correlation needs at least two samples.")
        return
    mat = log2_matrix(pivot_df, id_col, sample_cols, is_log2=False)
    corr = sample_correlation(mat)
    if corr.isna().all().all():
        st.info("No two samples share enough measured proteins to compare them.")
        return
    source = correlation_source(corr)
    groups = {s: sample_group_map[s] for s in corr.columns if s in sample_group_map}
    shared = shared_proteins(mat)
    head = heading(note=f"Each pair of samples is compared over the proteins measured in both: "
                        f"{median_shared_proteins(shared):,} for the median pair.")
    settings = library.settings(
        "proteins", sample_noun="Sample", metadata=sample_metadata(groups, list(corr.columns)),
        group_colors=group_colors(groups.values()),
        **({"title": head["title"]} if head.get("title") else {}),
    )
    heatmap = component(
        ClusteredHeatmap, _default_cache_dir(), cache_name("quantms_sample_correlation", source),
        lambda: {"data": source.lazy()},
        **only_accepted(ClusteredHeatmap, row_cluster=True, col_cluster=True, **settings),
    )
    height = max(500, min(1200, 250 + 28 * len(corr)))
    heatmap(state_manager=shared_state(), height=height, key=scoped("sample_correlation"), help=head.get("help"))
    medians = median_correlation(corr).dropna()
    outliers = least_similar(corr)
    if outliers:
        st.warning(
            "Less alike than the other samples: "
            + ", ".join(f"**{s}** (median r {medians[s]:.2f})" for s in outliers)
            + f", against a median of {medians.median():.2f}. Check whether these runs worked and belong to their group."
        )
    elif len(medians):
        st.caption(f"Median correlation between samples: {medians.median():.2f}, "
                   f"lowest {medians.iloc[0]:.2f} ({medians.index[0]}).")
