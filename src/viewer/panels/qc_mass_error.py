"""Precursor mass error of the target PSMs."""
from src.viewer.panels import qc


def render(ctx) -> None:
    qc.mass_error(ctx)
