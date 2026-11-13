"""Distribution of the p-values of the group comparison, as a Viewer component."""

from src.viewer.comparison_plots import comparison, draw_pvalues


def render(ctx) -> None:
    got = comparison(ctx)
    if got is None:
        return
    steps, _, _ = got
    draw_pvalues(steps["statistics_df"], "comparison-p-values")
