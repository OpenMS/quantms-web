"""The precursor XIC of the selected PSM, drawn with an OpenMS-Insight SeriesPlot."""
from __future__ import annotations

import shutil
from pathlib import Path

import polars as pl
import streamlit as st
from openms_insight import SeriesPlot

from src.common import precursor_xic as xic
from src.common.spectrum_view import psm_lookup

KEEP_CACHES = 30  # XIC plots kept on disk per run; older ones are removed
ISOTOPE_COLORS = {"M": "#1f77b4", "M+1": "#ff7f0e", "M+2": "#2ca02c", xic.MS2_SERIES: "#d62728"}


@st.cache_data(show_spinner="Tracing the precursor in the MS1 spectra...", max_entries=200)
def _isotope_traces(mzml: str, mtime: float, mz: float, charge: int, rt: float) -> object:
    # mtime only keys the cache
    return xic.extract_xic(mzml, mz, charge, rt)


@st.cache_data(show_spinner="Tracing the precursor in every run...", max_entries=50)
def _run_traces(mzmls: tuple, mtimes: tuple, mz: float, charge: int, rt: float) -> object:
    return xic.xic_across_runs({Path(p).stem: p for p in mzmls}, mz, charge, rt)


def _plot(cache_dir: Path, cache_id: str, data, series_colors: dict, title: str):
    """The SeriesPlot of ``data``, cached under ``cache_dir/xic``; the oldest caches beyond KEEP_CACHES go."""
    root = cache_dir / "xic"
    root.mkdir(parents=True, exist_ok=True)
    if not (root / cache_id).is_dir():
        caches = sorted((p for p in root.iterdir() if p.is_dir()), key=lambda p: p.stat().st_mtime)
        for stale in caches[: max(0, len(caches) - KEEP_CACHES + 1)]:
            shutil.rmtree(stale, ignore_errors=True)
        return SeriesPlot(
            cache_id=cache_id, data=pl.from_pandas(data).lazy(), cache_path=str(root),
            x_column="rt", y_column="intensity", series_column="series",
            series_names={name: name for name in series_colors},
            series_colors=series_colors, mode="lines+markers", title=title,
            x_label="Retention time (s)", y_label="Intensity", y_range_mode="nonnegative",
            svg_filename=cache_id,
        )
    return SeriesPlot(cache_id=cache_id, cache_path=str(root))


def show_precursor_xic(cache_dir: Path, idxml: Path, state_manager, *, mzml_files: list, key: str) -> None:
    """XIC of the precursor of the PSM selected in a run's PSM table (state: identification, file).

    ``mzml_files`` are the run's input mzML files in the order of the viewer
    caches (``psm_mzml_files(params)``), the MS1 spectra being read from them.
    """
    psm_id = state_manager.get_selection("identification")
    file_index = state_manager.get_selection("file")
    hits = psm_lookup(str(idxml), idxml.stat().st_mtime)
    if psm_id is None or file_index is None or int(psm_id) not in hits.index:
        st.info("Select a PSM in the table to see how its precursor eluted.")
        return
    hit = hits.loc[int(psm_id)]
    file_index = int(file_index)
    mzml = mzml_files[file_index] if file_index < len(mzml_files) else None
    if mzml is None or not Path(mzml).is_file():
        st.info(
            "The input mzML of this PSM is not in the workspace (a run opened from a ZIP keeps only the identified "
            "MS2 spectra), so its precursor cannot be traced over the MS1 spectra."
        )
        return
    charge, ms2_rt = int(hit["charge"]), float(hit["rt"])
    mz = xic.precursor_mz(hit["sequence"], charge, fallback=float(hit["mz"]))
    present = [Path(p) for p in mzml_files if p is not None and Path(p).is_file()]
    every_run = False
    if len(present) > 1:
        every_run = st.segmented_control(
            "Show", ["Isotope peaks in this run", "Monoisotopic peak in every run"],
            default="Isotope peaks in this run", key=f"{key}_xic_mode",
        ) == "Monoisotopic peak in every run"
    if every_run:
        traces = _run_traces(tuple(str(p) for p in present), tuple(p.stat().st_mtime for p in present),
                             mz, charge, ms2_rt)
        palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
                   "#17becf", "#aec7e8", "#ffbb78", "#98df8a"]
        names = [p.stem for p in present]
        colors = {name: palette[i % len(palette)] for i, name in enumerate(names)}
        colors[xic.MS2_SERIES] = "#d62728"
        data = xic.with_ms2_marker(traces, ms2_rt)
        st.caption(
            f"m/z {mz:.4f} ({charge}+) within ±{xic.PPM:g} ppm, ±{xic.RUNS_RT_WINDOW:g} s around the MS2 spectrum "
            "(the runs' retention times are not aligned). Similar peaks in the runs of a group are what the "
            "quantification compares."
        )
        absent = [name for name in names if name not in set(traces["series"])]
        if absent:
            st.caption(f"No MS1 spectra in this window: {', '.join(absent)}.")
        cache_id = f"xic_runs_{idxml.stem}_{int(psm_id)}"
        title = f"{hit['sequence']} {charge}+ in every run"
    else:
        traces = _isotope_traces(str(mzml), Path(mzml).stat().st_mtime, mz, charge, ms2_rt)
        data = xic.with_ms2_marker(traces, ms2_rt)
        colors = {xic.isotope_label(k): ISOTOPE_COLORS[xic.isotope_label(k)] for k in range(xic.ISOTOPES)}
        colors[xic.MS2_SERIES] = ISOTOPE_COLORS[xic.MS2_SERIES]
        summary = xic.summarize(traces, ms2_rt)
        cols = st.columns(4)
        cols[0].metric("Precursor m/z", f"{mz:.4f}", help=f"Monoisotopic, {charge}+, from the peptide sequence.")
        cols[1].metric("Apex RT (s)", "–" if summary["apex_rt"] is None else f"{summary['apex_rt']:.1f}")
        cols[2].metric(
            "MS2 vs. apex (s)", "–" if summary["ms2_offset"] is None else f"{summary['ms2_offset']:+.1f}",
            help="When the MS2 spectrum was taken, relative to the top of the elution peak.",
        )
        corr = summary["isotope_correlation"]
        cols[3].metric(
            "M / M+1 correlation", "–" if corr is None else f"{corr:.2f}",
            help="How closely the first isotope peak follows the monoisotopic one. Close to 1 for a clean "
            "peak; lower values point to another ion at the same m/z.",
        )
        st.caption(
            f"Monoisotopic and isotope peaks within ±{xic.PPM:g} ppm, ±{xic.RT_WINDOW:g} s around the MS2 "
            "spectrum (red line). The traces should rise and fall together."
        )
        cache_id = f"xic_{idxml.stem}_{int(psm_id)}"
        title = f"{hit['sequence']} {charge}+"
    if traces.empty:
        st.info("No MS1 spectra in this retention time window.")
        return
    _plot(cache_dir, cache_id, data, colors, title)(key=f"{key}_xic_plot", state_manager=state_manager, height=380)
