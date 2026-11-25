"""Single run: the PSMs of the run after FDR, down to the annotated spectrum of each."""
from src.viewer.panels import filtered


def render(ctx) -> None:
    filtered.render(ctx)
