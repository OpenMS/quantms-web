"""Where missing values occur, as a Viewer component."""

from src.common.postprocessing_plots import missingness_vs_intensity
from src.viewer.comparison_plots import comparison
from src.viewer.scope import show_fig


def render(ctx) -> None:
    got = comparison(ctx)
    if got is None:
        return
    steps, samples, _ = got
    show_fig(missingness_vs_intensity(steps["filtered_df"], samples), "comparison-missing-values")
