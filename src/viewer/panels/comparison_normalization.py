"""Sample intensities before and after normalization, as a Viewer component."""

from src.viewer.comparison_plots import comparison, draw_distributions


def render(ctx) -> None:
    got = comparison(ctx)
    if got is None:
        return
    steps, samples, _ = got
    draw_distributions(
        steps["imputed_df"], steps["normalized_df"], samples,
        after_is_log=steps.get("normalized_transform", "log2") in ("log2", "log10"),
        name="comparison-normalization",
    )
