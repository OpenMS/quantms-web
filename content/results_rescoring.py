"""Rescoring (Percolator) Results Page."""
import streamlit as st
from pathlib import Path
from src.common.common import page_setup
from src.common.results_helpers import get_workflow_dir
from openms_insight import Table, Heatmap, LinePlot, SequenceView, StateManager

params = page_setup()
st.title("Rescoring Results")

st.markdown(
    """
View PSMs after **Percolator** statistical validation. Percolator uses machine learning
to re-score PSMs and estimate false discovery rates (FDR) for more accurate results.
Click on a PSM to view the annotated spectrum and peptide sequence.
"""
)

st.info(
    "**Score:** The Posterior Error Probability (PEP) represents the probability that this PSM "
    "is incorrect. Lower values indicate higher confidence identifications."
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

workflow_dir = get_workflow_dir(st.session_state["workspace"])
perc_dir = workflow_dir / "results" / "percolator_results"
cache_dir = workflow_dir / "results" / "insight_cache"

if not perc_dir.exists():
    st.info("No rescoring results available yet. Please run the workflow first.")
    st.page_link("content/workflow_run.py", label="Go to Run", icon="🚀")
    st.stop()

perc_files = sorted(perc_dir.glob("*.idXML"))

if not perc_files:
    st.warning("No rescoring output files found.")
    st.stop()

selected_file = st.selectbox(
    "Select rescoring result file",
    perc_files,
    format_func=lambda x: x.name
)

cache_id_prefix = selected_file.stem

# Check if cache exists
if not (cache_dir / f"table_{cache_id_prefix}").is_dir():
    st.warning("Visualization cache not found. Please re-run the workflow.")
    st.stop()

# Initialize state manager for cross-component linking
state_manager = StateManager()

# Load components from cache (no data parameter needed)
table = Table(cache_id=f"table_{cache_id_prefix}", cache_path=str(cache_dir))
heatmap = Heatmap(cache_id=f"heatmap_{cache_id_prefix}", cache_path=str(cache_dir))
seq_view = SequenceView(cache_id=f"seqview_{cache_id_prefix}", cache_path=str(cache_dir))
line_plot = LinePlot(cache_id=f"lineplot_{cache_id_prefix}", cache_path=str(cache_dir))

# Render components
st.subheader("PSM Overview")
heatmap(state_manager=state_manager, height=350)

st.subheader("PSM Table")
table(state_manager=state_manager, height=533)

st.subheader("Peptide Sequence")
seq_view(key=f"seqview_{cache_id_prefix}", state_manager=state_manager, height=533)

st.subheader("MS2 Spectrum")
line_plot(key=f"lineplot_{cache_id_prefix}", state_manager=state_manager, height=450, sequence_view_key=f"seqview_{cache_id_prefix}")

st.markdown("---")
st.markdown("**Next step:** View filtered PSMs")
st.page_link("content/results_filtered.py", label="Go to Filtered PSMs", icon="🎯")
