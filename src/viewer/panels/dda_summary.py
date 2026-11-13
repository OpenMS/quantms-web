"""Identification summary: files, PSMs, peptides, proteins, ID rate and mass error as tiles."""
from src.viewer.panels import overview


def render(ctx) -> None:
    overview.dda_tiles(ctx)
