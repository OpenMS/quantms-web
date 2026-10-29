"""Identification QC Page."""
from pathlib import Path

import pandas as pd
import streamlit as st

from src.common.common import page_setup, show_fig
from src.common.qc_plots import psm_qc_table, precursor_error_histogram, psms_per_file, score_distribution
from src.common.results_helpers import get_workflow_dir
from src.workflow.ParameterManager import ParameterManager

params = page_setup()

st.title("Identification QC")
st.markdown(
    """
Quality checks of the identifications across all runs: how well targets separate
from decoys, how well the precursor masses are calibrated, and whether one run
identified clearly fewer spectra than the others.
"""
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

workflow_dir = get_workflow_dir(st.session_state["workspace"])
results_dir = workflow_dir / "results"
workflow_params = ParameterManager(workflow_dir, "TOPP Workflow").get_parameters_from_json()
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
    return load_qc_table(
        tuple(str(f) for f in files), tuple(f.stat().st_mtime for f in files), tuple(mzml_names)
    )


# All PSMs before filtering, so decoys and the FDR cutoff are visible
all_psms = qc_table(results_dir / "percolator_results")
if all_psms is None:
    all_psms = qc_table(results_dir / "comet_results")
filtered = qc_table(results_dir / "psm_filter")

if all_psms is None and filtered is None:
    st.info("No identification results available yet. Please run the workflow first.")
    st.page_link("content/workflow_run.py", label="Go to Run", icon="🚀")
    st.stop()

if all_psms is not None and not all_psms.empty:
    st.subheader("Score distribution")
    st.caption(
        "All PSMs before FDR filtering. Decoys model the incorrect matches: targets to the "
        "right of the decoy peak are the correct identifications. The dashed line is the "
        "lowest target score that passes the PSM FDR."
    )
    if not all_psms["decoy"].any():
        st.warning("No decoy PSMs found, so the score distribution cannot show the incorrect matches.")
    show_fig(score_distribution(all_psms, fdr), "qc-score-distribution")

if filtered is not None and not filtered.empty:
    st.subheader("Precursor mass error")
    st.caption(
        f"Target PSMs at the PSM FDR ({fdr * 100:.3g} %). A median far from 0 ppm points to a "
        "calibration offset; a distribution as wide as the precursor tolerance suggests the "
        "tolerance could be tighter."
    )
    show_fig(precursor_error_histogram(filtered), "qc-precursor-error")

    st.subheader("PSMs per file")
    st.caption(f"Target PSMs at the PSM FDR ({fdr * 100:.3g} %).")
    show_fig(psms_per_file(filtered), "qc-psms-per-file")
elif all_psms is not None:
    st.info("No filtered PSMs yet: the precursor error and PSMs per file appear once filtering has run.")
