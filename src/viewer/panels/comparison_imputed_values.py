"""The values filled in for missing ones, as a Viewer component."""

import streamlit as st

from src.common.postprocessing_plots import imputed_value_preview
from src.viewer.comparison_plots import comparison
from src.viewer.scope import show_fig


def render(ctx) -> None:
    got = comparison(ctx)
    if got is None:
        return
    steps, samples, _ = got
    st.caption(
        "Filled-in values (red) against observed ones (blue). They should sit at the low edge of the "
        "observed values; a separate spike far below them inflates fold changes."
    )
    show_fig(
        imputed_value_preview(steps["filtered_df"], steps["imputed_df"], samples), "comparison-imputed-values",
    )
