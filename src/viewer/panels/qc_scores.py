"""Target and decoy PSMs per score."""
from src.viewer.panels import qc


def render(ctx) -> None:
    qc.scores(ctx)
