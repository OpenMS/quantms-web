"""Distribution of the p-values of the group comparison, as a Viewer component."""

import streamlit as st

from src.common.postprocessing_plots import pvalue_histogram
from src.viewer.comparison_plots import comparison
from src.viewer.scope import show_fig


def render(ctx) -> None:
    got = comparison(ctx)
    if got is None:
        return
    steps, _, _ = got
    st.caption(
        "Raw p-values of all proteins. Flat with a peak near 0: the test fits and some proteins change. "
        "A peak near 1 or a U-shape points to the wrong test or non-log input."
    )
    show_fig(pvalue_histogram(steps["statistics_df"]), "comparison-p-values")
