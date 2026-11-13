"""The values filled in for missing ones, as a Viewer component."""

from src.viewer.comparison_plots import comparison, draw_imputed


def render(ctx) -> None:
    got = comparison(ctx)
    if got is None:
        return
    steps, samples, _ = got
    draw_imputed(steps["filtered_df"], steps["imputed_df"], samples, "comparison-imputed-values")
