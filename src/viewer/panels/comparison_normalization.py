"""Sample intensities before and after normalization, as a Viewer component."""

import streamlit as st

from src.common.postprocessing_plots import sample_distributions
from src.viewer.comparison_plots import comparison
from src.viewer.scope import show_fig


def render(ctx) -> None:
    got = comparison(ctx)
    if got is None:
        return
    steps, samples, sample_group_map = got
    st.caption(
        "Each sample's intensities before (top) and after (bottom) normalization. After a good "
        "normalization the medians line up."
    )
    show_fig(
        sample_distributions(
            steps["imputed_df"], steps["normalized_df"], samples, sample_group_map,
            after_is_log=steps.get("normalized_transform", "log2") in ("log2", "log10"),
        ),
        "comparison-normalization",
    )
