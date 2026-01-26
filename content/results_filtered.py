"""Filtered PSMs Results Page."""
import streamlit as st
from pathlib import Path
from src.common.common import page_setup
from src.common.results_helpers import get_workflow_dir
from openms_insight import Table, Heatmap, LinePlot, SequenceView, StateManager

params = page_setup()
st.title("Filtered PSMs")

st.markdown(
    """
View FDR-controlled peptide identifications after **IDFilter** processing.
These are high-confidence PSMs that passed the specified FDR threshold.
Click on a PSM to view the annotated spectrum and peptide sequence.
"""
)

st.info("Explore the PSM scatterplot along with the detailed PSM table.")

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

workflow_dir = get_workflow_dir(st.session_state["workspace"])
filter_dir = workflow_dir / "results" / "filter_results"
cache_dir = workflow_dir / "results" / "insight_cache"

if not filter_dir.exists():
    st.info("No filtered results available yet. Please run the workflow first.")
    st.page_link("content/workflow_run.py", label="Go to Run", icon="🚀")
    st.stop()

filter_files = sorted(filter_dir.glob("*.idXML"))

if not filter_files:
    st.warning("No filtering output files found.")
    st.stop()

selected_file = st.selectbox(
    "Select filtering result file",
    filter_files,
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
table(state_manager=state_manager, height=400)

st.subheader("Peptide Sequence")
seq_view(key=f"seqview_{cache_id_prefix}", state_manager=state_manager, height=400)

st.subheader("MS2 Spectrum")
line_plot(key=f"lineplot_{cache_id_prefix}", state_manager=state_manager, height=450, sequence_view_key=f"seqview_{cache_id_prefix}")

st.markdown("---")
st.markdown("**Next step:** View abundance quantification")
st.page_link("content/results_abundance.py", label="Go to Abundance", icon="📋")
