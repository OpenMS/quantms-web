"""Download Results Page."""
import streamlit as st
from src.common.common import page_setup
from src.common.guide import page_intro
from src.WorkflowTest import WorkflowTest

params = page_setup()

wf = WorkflowTest()

st.title("Download results")

page_intro(
    """
Every file the analysis produced, singly or as one ZIP. The ones most people need:

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
