"""The MS1 peak map of a file with the precursors of its PSMs, and its chromatograms."""

import hashlib
import shutil
from pathlib import Path

import polars as pl
import streamlit as st
from openms_insight import Heatmap, SeriesPlot

from src.common import peak_map
from src.common.results_helpers import extract_filename_from_idxml, get_run_dir, get_run_params
from src.common.spectrum_view import psm_lookup, psm_mzml_files
from src.viewer.scope import PSM_PRIVATE, needs_results, scoped, shared_state, zoom

# Owner of the PSM row ids shared with the PSM table (src/viewer/panels/filtered.py)
PSM_OWNER = "filtered"
CHROMATOGRAM_COLORS = {peak_map.TIC: "#3157e9", peak_map.BPC: "#ff7501"}


@st.cache_data(show_spinner="Reading the MS1 spectra...", max_entries=4)
def _peaks(mzml: str, mtime: float, idxml: str, id_mtime: float) -> tuple:
    """MS1 peaks with the PSM precursors marked, and the chromatograms (mtimes only key the cache)."""
    peaks = peak_map.ms1_peaks(mzml)
    psms = pl.from_pandas(psm_lookup(idxml, id_mtime).reset_index()[["id_idx", "mz", "rt"]])
    marked = peak_map.mark_precursors(peaks, psms)
    return marked, peak_map.chromatograms(peaks)


def _cache_id(prefix: str, *paths: Path) -> str:
    stamp = "".join(f"{p}{p.stat().st_mtime}" for p in paths)
    return f"{prefix}_{hashlib.sha1(stamp.encode()).hexdigest()[:10]}"


def _synced_file(key: str, idxmls: list, names: list, state_manager) -> None:
    """Show the file of a PSM picked in another component, once per pick."""
    file_index = state_manager.get_selection("file")
    synced = f"{key}_synced"
    if file_index is None or st.session_state.get(synced) == file_index:
        return
    st.session_state[synced] = file_index
    if int(file_index) >= len(names):
        return
    for idxml in idxmls:
        if extract_filename_from_idxml(idxml, names) == names[int(file_index)]:
            st.session_state[key] = idxml
            return


def render(ctx) -> None:
    run_dir = get_run_dir(ctx.workspace)
    run_params = get_run_params(ctx.workspace)
    filter_dir = run_dir / "results" / "psm_filter"
    cache_root = run_dir / "results" / "insight_cache" / "peak_map"
    idxmls = sorted(filter_dir.glob("*.idXML")) if filter_dir.is_dir() else []
    if not idxmls:
        return needs_results()
    mzmls = psm_mzml_files(run_params)
    names = [Path(f).name for f in mzmls]
    state_manager = shared_state(PSM_OWNER, PSM_PRIVATE)
    key = scoped("peak_map_file")
    _synced_file(key, idxmls, names, state_manager)
    if st.session_state.get(key) not in idxmls:
        st.session_state.pop(key, None)
    idxml = st.selectbox(
        "File", idxmls, key=key,
        format_func=lambda p: Path(extract_filename_from_idxml(p, names)).stem,
    )
    source = extract_filename_from_idxml(idxml, names)
    mzml = next((Path(f) for f in mzmls if Path(f).name == source), None)
    if mzml is None or not mzml.is_file():
        st.info(
            "The input mzML of this file is not in the workspace (a run opened from a ZIP keeps only the "
            "identified MS2 spectra), so its MS1 spectra cannot be shown."
        )
        return

    cache_root.mkdir(parents=True, exist_ok=True)
    map_id = _cache_id(f"peakmap_{idxml.stem}", mzml, idxml)
    chrom_id = _cache_id(f"chromatogram_{idxml.stem}", mzml)
    if not (cache_root / map_id).is_dir() or not (cache_root / chrom_id).is_dir():
        marked, chromatogram = _peaks(str(mzml), mzml.stat().st_mtime, str(idxml), idxml.stat().st_mtime)
        for stale in cache_root.glob(f"peakmap_{idxml.stem}_*"):
            if stale.name != map_id:
                shutil.rmtree(stale, ignore_errors=True)
        heatmap = Heatmap(
            cache_id=map_id, data=marked.lazy(), cache_path=str(cache_root),
            x_column="rt", y_column="mz", intensity_column="intensity",
            interactivity={"identification": peak_map.PRECURSOR_COLUMN},
            highlight={"identification": peak_map.PRECURSOR_COLUMN},
            title=f"MS1 peak map: {mzml.stem}", x_label="Retention time (s)", y_label="m/z",
            intensity_label="Intensity", svg_filename=f"{mzml.stem}_peak_map",
            **zoom(Heatmap, "peak_map"),
        )
        series = SeriesPlot(
            cache_id=chrom_id, data=chromatogram.lazy(), cache_path=str(cache_root),
            x_column="rt", y_column="intensity", series_column="series",
            series_names={name: name for name in CHROMATOGRAM_COLORS}, series_colors=CHROMATOGRAM_COLORS,
            mode="lines", title=f"Chromatograms: {mzml.stem}", x_label="Retention time (s)",
            y_label="Intensity", y_range_mode="nonnegative", svg_filename=f"{mzml.stem}_chromatograms",
            **zoom(SeriesPlot, "peak_map_chromatogram"),
        )
    else:
        heatmap = Heatmap(cache_id=map_id, cache_path=str(cache_root), **zoom(Heatmap, "peak_map"))
        series = SeriesPlot(cache_id=chrom_id, cache_path=str(cache_root), **zoom(SeriesPlot, "peak_map_chromatogram"))

    st.caption(
        "Every MS1 peak of the file by retention time and m/z. A click on the peak a PSM's MS2 spectrum was "
        "taken from picks that PSM; the PSM picked in the table is ringed."
    )
    heatmap(key=scoped("peak_map_plot"), state_manager=state_manager, height=450)
    series(key=scoped("peak_map_chromatogram_plot"), state_manager=state_manager, height=300)
