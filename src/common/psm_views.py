"""The viewer caches of one identification file: PSM table, plot, sequence and spectrum.

The workflow builds them as it runs. A run opened from a ZIP has no caches
(the ZIP leaves them out); ``rebuild_psm_views`` makes them again from the
run's own spectra subset, or from the input mzML when it is there.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl
from openms_insight import Heatmap, LinePlot, SequenceView, Table

from src.common.identified_spectra import has_spectra, spectra_sources
from src.common.results_helpers import (
    build_spectra_cache, extract_filename_from_idxml, parse_idxml, peaks_for_identifications, psm_table_columns,
    write_psm_exports,
)

CACHE_DIR = "insight_cache"
# Folder of the raw database search hits; their score (ProSE's ln(hyperscore))
# is better when higher, unlike Percolator's and IDFilter's (and earlier runs'
# Comet e-value in comet_results/), which are better when lower
SEARCH_RESULTS = "search_results"


def score_sort(idxml_path: Path) -> str:
    """The PSM table's initial sort on ``score`` for this idXML, best matches first."""
    return "desc" if Path(idxml_path).parent.name == SEARCH_RESULTS else "asc"


def build_psm_views(idxml_path: Path, cache_dir: Path, spectra_df: pl.DataFrame | dict, filename_to_index: dict,
                    tolerance: float, ppm: bool, export_dir: Path | None = None) -> bool:
    """Cache the Table, Heatmap, SequenceView and LinePlot of one idXML.

    Args:
        spectra_df: MS2 peaks, or the per-file cache of ``write_spectra_cache``.
        export_dir: Also write the PSM and fragment ion TSV exports there.

    Returns:
        False when the idXML holds no identification (nothing is cached).
    """
    idxml_path = Path(idxml_path)
    prefix = idxml_path.stem
    id_df, _ = parse_idxml(idxml_path, filename_to_index)
    if id_df.is_empty():
        return False

    Table(
        cache_id=f"table_{prefix}",
        data=id_df.lazy(),
        cache_path=str(cache_dir),
        interactivity={"file": "file_index", "spectrum": "scan_id", "identification": "id_idx"},
        column_definitions=psm_table_columns(),
        initial_sort=[{"column": "score", "dir": score_sort(idxml_path)}],
        index_field="id_idx",
    )
    Heatmap(
        cache_id=f"heatmap_{prefix}",
        data=id_df.lazy(),
        cache_path=str(cache_dir),
        x_column="rt",
        y_column="mz",
        intensity_column="score",
        interactivity={"identification": "id_idx"},
    )
    seq_view = SequenceView(
        cache_id=f"seqview_{prefix}",
        sequence_data=id_df.lazy().select(["id_idx", "sequence", "charge", "file_index", "scan_id"]).rename({
            "id_idx": "sequence_id",
            "charge": "precursor_charge",
        }),
        peaks_data=peaks_for_identifications(spectra_df, id_df),
        filters={
            "identification": "sequence_id",
            "file": "file_index",
            "spectrum": "scan_id",
        },
        interactivity={"peak": "peak_id"},
        cache_path=str(cache_dir),
        deconvolved=False,
        annotation_config={
            "ion_types": ["b", "y"],
            "neutral_losses": True,
            "tolerance": tolerance,
            "tolerance_ppm": ppm,
        },
    )
    LinePlot.from_sequence_view(
        seq_view,
        cache_id=f"lineplot_{prefix}",
        cache_path=str(cache_dir),
        title="Annotated Spectrum",
        styling={
            "unhighlightedColor": "#CCCCCC",
            "highlightColor": "#E74C3C",
            "selectedColor": "#F3A712",
        },
    )
    if export_dir is not None:
        # Downloadable exports: every PSM, and every PSM's matched fragment
        # ions against its own spectrum
        write_psm_exports(id_df, seq_view, Path(export_dir), prefix)
    return True


def rebuild_psm_views(results_dir: Path, idxml_path: Path, mzml_files: list, tolerance: float, ppm: bool) -> bool:
    """Make the viewer caches of ``idxml_path`` again when they are missing.

    Only runs that carry their identified spectra (``results/spectra``) are
    rebuilt; for others nothing changes. Spectra come from the input mzML
    when it is there, else from the run's subset of it.

    Args:
        mzml_files: The run's input mzML paths, in the order that numbers them.

    Returns:
        Whether the caches exist afterwards.
    """
    cache_dir = Path(results_dir, CACHE_DIR)
    idxml_path = Path(idxml_path)
    if (cache_dir / f"table_{idxml_path.stem}").is_dir():
        return True
    if not mzml_files or not has_spectra(results_dir):
        return False
    names = [Path(f).name for f in mzml_files]
    filename_to_index = {name: i for i, name in enumerate(names)}
    # One idXML holds one input file; read only that file's spectra
    source_name = extract_filename_from_idxml(idxml_path, names)
    sources = dict(zip(names, spectra_sources(mzml_files, results_dir)))
    source = sources.get(source_name)
    if source is None:
        return False
    spectra_df, filename_to_index = build_spectra_cache(source.parent, filename_to_index, [source])
    cache_dir.mkdir(parents=True, exist_ok=True)
    return build_psm_views(idxml_path, cache_dir, spectra_df, filename_to_index, tolerance, ppm)
