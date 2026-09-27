"""Download Results Page."""
import streamlit as st
from src.common.common import page_setup
from src.WorkflowTest import WorkflowTest

params = page_setup()

wf = WorkflowTest()

st.title("Download Results")

st.markdown(
    """
Every file the workflow produced: identifications (idXML) from each search step,
quantification tables, the spectral library if one was built, and the exports in
**exports/**: all PSMs per run (`*_psms.tsv`) and every PSM's matched fragment
ions per scan (`*_fragment_ions.tsv`).
"""
)

wf.show_results_download_section(exclude=["insight_cache"])
