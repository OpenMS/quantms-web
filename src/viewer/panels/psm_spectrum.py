"""Annotated spectrum of the picked PSM."""
from src.viewer import protein_data
from src.viewer.scope import needs_results


def render(ctx) -> None:
    data = protein_data.load(ctx)
    if data is None:
        return needs_results()
    protein_data.psm_spectrum(data)
