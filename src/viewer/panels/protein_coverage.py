"""Sequence coverage of the picked protein: its peptides drawn where they lie
on its sequence, one track per condition; a click picks the peptide."""

from src.viewer import protein_data
from src.viewer.scope import needs_results


def render(ctx) -> None:
    data = protein_data.load(ctx)
    if data is None:
        return needs_results()
    protein_data.protein_coverage(data)
