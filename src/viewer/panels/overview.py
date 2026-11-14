"""results_overview as a Viewer component."""

from pathlib import Path
import streamlit as st
from src.common.empty_result import fdr_summary
from src.common.insight_tables import show_insight_table
from src.common.overview import PER_FILE_COLUMNS, run_overview
from src.common.results_helpers import get_run_dir, get_run_params, search_tolerances
from src.viewer.scope import needs_results


@st.cache_data(show_spinner="Summarizing the run...")
def _overview(workflow_dir: str, files: tuple, precursor_ppm, mtimes: tuple):
    """run_overview, cached on the inputs' and results' modification times."""
    return run_overview(workflow_dir, list(files), precursor_ppm=precursor_ppm)


def load_overview(ctx):
    """The shown run's per-file table and key figures, with what they were read from.

    Returns (per_file, figures, workflow_params, cache_dir, sources), or None
    when the run has no filtered identifications. Shared by the run summary and
    the diagnostics component, so the files are read once.
    """
    run_dir = get_run_dir(ctx.workspace)
    filter_dir = run_dir / "results" / "psm_filter"
    workflow_params = get_run_params(ctx.workspace)
    mzml_files = [str(f) for f in workflow_params.get("mzML-files", [])]
    if not filter_dir.is_dir() or not any(filter_dir.glob("*.idXML")):
        return None
    sources = sorted(filter_dir.glob("*.idXML")) + [run_dir / "logs" / "commands-and-run-times.log"]
    tolerances = search_tolerances(workflow_params)
    precursor_ppm = tolerances["precursor"] if tolerances["precursor_unit"] == "ppm" else None
    mtimes = tuple(Path(p).stat().st_mtime if Path(p).exists() else 0.0 for p in mzml_files + sources)
    per_file, figures = _overview(str(run_dir), tuple(mzml_files), precursor_ppm, mtimes)
    return per_file, figures, workflow_params, run_dir / "results" / "insight_cache", sources


def dda_tiles(ctx) -> None:
    """Identification summary as tiles, with the run's filters as a grey line."""
    from src.viewer import protein_data

    loaded = load_overview(ctx)
    if loaded is None:
        return needs_results()
    _, figures, workflow_params, _, _ = loaded
    data = protein_data.load(ctx)
    proteins = protein_data.protein_count(data) if data is not None else None
    fdr = float(workflow_params.get("psm-fdr-percent", 1.0))
    tiles = [
        ("Files", _show(figures["Files"]), f"Filters of this run: {fdr_summary(workflow_params)}."),
        ("Target PSMs", _show(figures["Target PSMs"]), f"At {fdr:g} % PSM FDR, over all files."),
        ("Peptides", _show(figures["Peptides"]), "Distinct target peptide sequences."),
        ("Proteins", _show(proteins), "Identified target proteins, without contaminants."),
        ("PSMs per MS2 spectrum", _show(figures["PSMs per MS2 (%)"], " %"), None),
        ("Median precursor error", _show(figures["Median precursor error (ppm)"], " ppm"), None),
    ]
    for col, (label, value, help_text) in zip(st.columns(len(tiles)), tiles):
        col.metric(label, value, help=help_text)


def files_table(ctx) -> None:
    """One row per file: spectra, PSMs, ID rate, peptides, mass error."""
    loaded = load_overview(ctx)
    if loaded is None:
        return needs_results()
    per_file, _, workflow_params, cache_dir, sources = loaded
    fdr = float(workflow_params.get("psm-fdr-percent", 1.0))
    # The fractional columns go as text, since the table stores floats as float32
    # and would show 36.8 as 36.79999923706055
    per_file_view = per_file[PER_FILE_COLUMNS].copy()
    for col in ("ID rate (%)", "Mass error (ppm)"):
        per_file_view[col] = per_file_view[col].map(lambda v: "" if v is None or v != v else f"{v:g}")
    show_insight_table(
        per_file_view,
        name="overview_files",
        cache_dir=cache_dir,
        source_files=sources,
        title="Files",
        note=f"Counts at {fdr:g} % PSM FDR.",
        height=300,
        # Share the page width between the columns rather than folding the last ones away
        layout="fitColumns",
    )


def _show(value, unit: str = "") -> str:
    return "–" if value is None else f"{value:,}{unit}"


def render(ctx) -> None:
    """The former Run summary (kept for saved views): the tiles and the files."""
    dda_tiles(ctx)
    files_table(ctx)
