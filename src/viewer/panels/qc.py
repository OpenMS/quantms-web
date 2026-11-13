"""Identification quality: score distribution, precursor mass error and PSMs per file."""

from pathlib import Path
import pandas as pd
import polars as pl
import streamlit as st
from openms_insight import SeriesPlot
from src.common.insight_tables import plot_cache
from src.common.qc_plots import (
    DECOY_COLOR, TARGET_COLOR, fdr_cutoff, mass_error_histogram, psm_qc_table, psms_per_file, score_histogram,
)
from src.common.results_helpers import get_run_dir, get_run_params, search_results_dir
from src.viewer.scope import needs_results, scoped, shared_state, show_fig, zoom


@st.cache_data(show_spinner="Reading identifications...")
def _load_qc_table(paths: tuple[str, ...], mtimes: tuple[float, ...], mzml_names: tuple[str, ...]) -> pd.DataFrame:
    # mtimes only key the cache, so a re-run of the workflow invalidates it
    return pd.concat([psm_qc_table(Path(p), list(mzml_names)) for p in paths], ignore_index=True)


def load_qc(ctx):
    """(all PSMs before filtering, PSMs after the FDR filter, PSM FDR as a fraction); None without results."""
    results_dir = get_run_dir(ctx.workspace) / "results"
    workflow_params = get_run_params(ctx.workspace)
    fdr = min(float(workflow_params.get("psm-fdr-percent", 1.0)), 100.0) / 100.0
    mzml_names = tuple(Path(f).name for f in workflow_params.get("mzML-files", []) or [])

    def table(directory: Path) -> pd.DataFrame | None:
        files = sorted(directory.glob("*.idXML"))
        if not files:
            return None
        return _load_qc_table(tuple(str(f) for f in files), tuple(f.stat().st_mtime for f in files), mzml_names)

    # All PSMs before filtering, so decoys and the FDR cutoff are visible
    all_psms = table(results_dir / "percolator_results")
    if all_psms is None:
        all_psms = table(search_results_dir(results_dir))
    filtered = table(results_dir / "psm_filter")
    if all_psms is None and filtered is None:
        return None
    return all_psms, filtered, fdr


def _series(name: str, data: pd.DataFrame, **kwargs) -> None:
    SeriesPlot(
        **plot_cache(f"quantms_{name}", data), **zoom(SeriesPlot, name),
        data=pl.from_pandas(data).lazy(), x_column="x", y_column="y", mode="lines", **kwargs,
    )(key=scoped(name), state_manager=shared_state(), height=340)


def scores(ctx) -> None:
    """Target and decoy PSMs per score, before the FDR filter."""
    loaded = load_qc(ctx)
    if loaded is None or loaded[0] is None or loaded[0].empty:
        return needs_results()
    all_psms, _, fdr = loaded
    data, label, _ = score_histogram(all_psms)
    cutoff = fdr_cutoff(all_psms, fdr)
    targets, decoys = int((~all_psms["decoy"]).sum()), int(all_psms["decoy"].sum())
    line = f"{targets:,} target and {decoys:,} decoy PSMs"
    if cutoff is not None:
        line += f" · PSM FDR {fdr * 100:.3g} % passes {label} ≥ {cutoff:.3g}"
    st.caption(line + ".")
    if not decoys:
        st.warning("No decoy PSMs found, so the score distribution cannot show the incorrect matches.")
    _series(
        "qc_scores", data, series_column="kind",
        series_names={"target": "Targets", "decoy": "Decoys"},
        series_colors={"target": TARGET_COLOR, "decoy": DECOY_COLOR},
        x_label=label, y_label="PSMs",
    )


def mass_error(ctx) -> None:
    """Precursor mass error of the target PSMs after the FDR filter."""
    loaded = load_qc(ctx)
    if loaded is None or loaded[1] is None or loaded[1].empty:
        return needs_results()
    _, filtered, fdr = loaded
    data, median = mass_error_histogram(filtered)
    if median is None:
        return needs_results("No precursor errors in this run.")
    st.caption(f"Target PSMs at {fdr * 100:.3g} % PSM FDR · median {median:+.2f} ppm, isotope errors removed.")
    _series(
        "qc_mass_error", data, series_colors=None, x_label="Precursor error (ppm)", y_label="PSMs",
        show_legend=False, y_range_mode="nonnegative",
    )


def render(ctx) -> None:
    """The former Identification quality component (kept for saved views)."""
    loaded = load_qc(ctx)
    if loaded is None:
        return needs_results()
    scores(ctx)
    mass_error(ctx)
    filtered = loaded[1]
    if filtered is not None and not filtered.empty:
        show_fig(psms_per_file(filtered), "qc-psms-per-file")
