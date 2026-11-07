"""Where missing values occur, as a Viewer component."""

import streamlit as st

from src.common.postprocessing_plots import missingness_vs_intensity
from src.viewer.comparison_plots import comparison
from src.viewer.scope import show_fig


def render(ctx) -> None:
    got = comparison(ctx)
    if got is None:
        return
    steps, samples, _ = got
    st.caption(
        "Mean intensity of proteins with and without missing values. A clearly lower orange curve "
        "means values fall below the detection limit (MNAR); overlapping curves mean random dropout (MAR)."
    )
    show_fig(missingness_vs_intensity(steps["filtered_df"], samples), "comparison-missing-values")
