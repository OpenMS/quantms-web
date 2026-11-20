"""One row per file of the run: spectra, PSMs, ID rate, peptides and mass error."""
from src.viewer.panels import overview


def render(ctx) -> None:
    overview.files_table(ctx)
