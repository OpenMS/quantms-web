"""results_qc as a Viewer component."""

from pathlib import Path
import pandas as pd
import streamlit as st
from src.common.qc_plots import psm_qc_table, precursor_error_histogram, psms_per_file, score_distribution
from src.common.results_helpers import get_run_dir, get_run_params, search_results_dir
from src.viewer.scope import needs_results, show_fig


def render(ctx) -> None:
    run_dir = get_run_dir(ctx.workspace)
    results_dir = run_dir / "results"
    workflow_params = get_run_params(ctx.workspace)
    fdr = min(float(workflow_params.get("psm-fdr-percent", 1.0)), 100.0) / 100.0
    mzml_names = [Path(f).name for f in workflow_params.get("mzML-files", []) or []]

    @st.cache_data(show_spinner="Reading identifications...")
    def load_qc_table(paths: tuple[str, ...], mtimes: tuple[float, ...], mzml_names: tuple[str, ...]) -> pd.DataFrame:
        # mtimes only key the cache, so a re-run of the workflow invalidates it
        return pd.concat([psm_qc_table(Path(p), list(mzml_names)) for p in paths], ignore_index=True)

    def qc_table(directory: Path) -> pd.DataFrame | None:
        files = sorted(directory.glob("*.idXML"))
        if not files:
            return None
        return load_qc_table(tuple(str(f) for f in files), tuple(f.stat().st_mtime for f in files), tuple(mzml_names))

    # All PSMs before filtering, so decoys and the FDR cutoff are visible
    all_psms = qc_table(results_dir / "percolator_results")
    if all_psms is None:
        all_psms = qc_table(search_results_dir(results_dir))
    filtered = qc_table(results_dir / "psm_filter")
    if all_psms is None and filtered is None:
        return needs_results()
    if all_psms is not None and not all_psms.empty:
        st.subheader("Score distribution")
        if not all_psms["decoy"].any():
            st.warning("No decoy PSMs found, so the score distribution cannot show the incorrect matches.")
        show_fig(score_distribution(all_psms, fdr), "qc-score-distribution")
    if filtered is not None and not filtered.empty:
        st.subheader("Precursor mass error")
        st.caption(f"Target PSMs at the PSM FDR ({fdr * 100:.3g} %).")
        show_fig(precursor_error_histogram(filtered), "qc-precursor-error")

        st.subheader("PSMs per file")
        st.caption(f"Target PSMs at the PSM FDR ({fdr * 100:.3g} %).")
        show_fig(psms_per_file(filtered), "qc-psms-per-file")
    elif all_psms is not None:
        st.info("No filtered PSMs yet: the precursor error and PSMs per file appear once filtering has run.")
