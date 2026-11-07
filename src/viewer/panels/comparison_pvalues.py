"""Distribution of the p-values of the group comparison, as a Viewer component."""

from src.common.postprocessing_plots import pvalue_histogram
from src.viewer.comparison_plots import comparison
from src.viewer.scope import show_fig


def render(ctx) -> None:
    got = comparison(ctx)
    if got is None:
        return
    steps, _, _ = got
    show_fig(pvalue_histogram(steps["statistics_df"]), "comparison-p-values")
