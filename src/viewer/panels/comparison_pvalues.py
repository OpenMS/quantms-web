"""Distribution of the p-values of the group comparison, as a Viewer component."""

from src.viewer.comparison_plots import comparison, draw_pvalues, plot_height


def render(ctx) -> None:
    got = comparison(ctx)
    if got is None:
        return
    steps, samples, _ = got
    draw_pvalues(steps["statistics_df"], "comparison-p-values", height=plot_height(samples))
