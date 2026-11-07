"""The values filled in for missing ones, as a Viewer component."""

from src.common.postprocessing_plots import imputed_value_preview
from src.viewer.comparison_plots import comparison
from src.viewer.scope import show_fig


def render(ctx) -> None:
    got = comparison(ctx)
    if got is None:
        return
    steps, samples, _ = got
    show_fig(
        imputed_value_preview(steps["filtered_df"], steps["imputed_df"], samples), "comparison-imputed-values",
    )
