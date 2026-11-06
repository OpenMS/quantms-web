"""Export, below the Viewer: the runs, the files of one run, its provenance and spectral library."""
from pathlib import Path
from types import SimpleNamespace

import streamlit as st

from src.viewer.panels import download, library


def export_section(workspace: Path) -> None:
    local = st.session_state.get("location") == "local"
    st.header("Export")
    st.markdown(
        "Every analysis you start is kept as its own run. 👁️ shows a run on Explore results, ✏️ renames it, "
        "📦 packs its results, settings, logs and identified spectra into one ZIP, and 🗑️ deletes it. "
        "Below the runs are the files of one run, singly or as one ZIP"
        + (", or exported to a folder." if local else ".")
    )
    with st.expander("📖 Which files most people need"):
        st.markdown(
            """
- **quant_results/**: the protein and peptide abundance table (CSV, opens in Excel or R).
- **exports/**: all peptide-spectrum matches per run (`*_psms.tsv`) and each match's
  fragment ions (`*_fragment_ions.tsv`).
- **psm_filter/**, **percolator_results/**, **search_results/**: identifications after each
  step in the open idXML format, for OpenMS and other proteomics software.
- **qpx/**: the peptide-spectrum matches in the QPX format (`run.psm.parquet`) with the run's
  provenance (`run.provenance.parquet`): every tool, its version and every parameter.
- The spectral library, if you built one.

The group comparison (fold changes and p-values) is downloaded with the download button of
the Group comparison table.
"""
        )
    ctx = SimpleNamespace(workspace=workspace, params={})
    download.render(ctx)
    with st.expander("Spectral library"):
        library.render(ctx)
