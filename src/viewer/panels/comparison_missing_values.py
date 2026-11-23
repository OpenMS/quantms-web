"""Where missing values occur, as a Viewer component."""

from src.viewer.comparison_plots import comparison, draw_missingness


def render(ctx) -> None:
    got = comparison(ctx)
    if got is None:
        return
    steps, samples, _ = got
    draw_missingness(steps["filtered_df"], samples, "comparison-missing-values")
