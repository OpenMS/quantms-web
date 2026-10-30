"""Download Results Page."""
import streamlit as st
from src.common.common import page_setup
from src.common.guide import page_intro
from src.WorkflowTest import WorkflowTest

params = page_setup()

wf = WorkflowTest()

st.title("Runs and downloads")

page_intro(
    """
Every analysis you start is kept as its own run, so you can change a setting, run again
and compare. The list below shows all runs: 👁️ shows a run on the results pages (the
**Showing run** box at the top of each results page does the same), ✏️ renames it, 📦 packs its results,
settings and logs into one ZIP to download, and 🗑️ deletes it.

Below the list are the files of one run, singly or as one ZIP. The ones most people need:

- **quant_results/**: the protein and peptide abundance table (CSV, opens in Excel or R).
- **exports/**: all peptide-spectrum matches per run (`*_psms.tsv`) and each match's
  fragment ions (`*_fragment_ions.tsv`).
- **psm_filter/**, **percolator_results/**, **comet_results/**: identifications after each
  step in the open idXML format, for OpenMS and other proteomics software.
- The spectral library, if you built one.

The group comparison (fold changes and p-values) is downloaded from the table on the
**Compare groups** page with the table's download button.
""",
)

wf.show_results_download_section(exclude=["insight_cache"])
