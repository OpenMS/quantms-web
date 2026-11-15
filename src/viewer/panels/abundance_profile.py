"""Abundance of the picked protein in every sample, by group (the template's
library panel); with several proteins picked together, each of them."""
import pandas as pd
import streamlit as st
from openms_insight import SeriesPlot

from src.common.differential import lead_accession
from src.common.insight_tables import _default_cache_dir
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from src.viewer.scope import GO_SET, heading, needs_results, scoped, shared_state, zoom
from src.viz import profile


def profile_source(pivot_df: pd.DataFrame, id_col: str, sample_group_map: dict) -> pd.DataFrame:
    """The library's contract: one row per protein (lead accession) and sample."""
    samples = [c for c in pivot_df.columns if c in sample_group_map]
    long = pivot_df[[id_col, *samples]].melt(id_vars=id_col, var_name="sample", value_name="value")
    return pd.DataFrame({
        "accession": long[id_col].map(lead_accession),
        "sample": long["sample"].astype(str),
        "group": long["sample"].map(sample_group_map).astype(str),
        "value": pd.to_numeric(long["value"], errors="coerce"),
    })


def profile_table(pivot_df: pd.DataFrame, id_col: str, sample_group_map: dict) -> pd.DataFrame:
    """The points and group means the panel draws (log2 abundance)."""
    samples = [c for c in pivot_df.columns if c in sample_group_map]
    return profile.profile_table(profile_source(pivot_df, id_col, sample_group_map), "accession",
                                 log2=True, samples=samples)


def render(ctx) -> None:
    result = get_abundance_data(ctx.workspace)
    if result is None:
        return needs_results()
    pivot_df, _, group_map = result
    id_col = get_id_column(ctx.workspace, pivot_df)
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    source = profile_source(pivot_df, id_col, sample_group_map)
    if not (source["value"] > 0).any():
        st.info("No abundances in this run.")
        return
    state_manager = shared_state()
    # Several proteins picked together (the volcano's box or lasso, the check
    # boxes of the group comparison): each of them as a series
    several = bool(profile.picked_set(state_manager, GO_SET))
    head = heading()
    profile.show_profile(
        source, "accession", GO_SET if several else "protein",
        "quantms_abundance_profiles" if several else "quantms_abundance_profile_v3",
        cache_path=_default_cache_dir(), key=scoped("abundance_profile"), state_manager=state_manager,
        zoom=zoom(SeriesPlot, "abundance_profile").get("zoom_identifier"), help=head.get("help"),
        title=head.get("title"), several=several,
        table_args={"log2": True, "samples": [c for c in pivot_df.columns if c in sample_group_map]},
        x_label="Sample", y_label="log2 abundance",
    )
