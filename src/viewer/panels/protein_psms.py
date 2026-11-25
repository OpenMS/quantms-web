"""PSMs of the picked protein and peptide."""
from src.viewer import protein_data
from src.viewer.scope import needs_results


def render(ctx) -> None:
    data = protein_data.load(ctx)
    if data is None:
        return needs_results()
    protein_data.protein_psms(data)
