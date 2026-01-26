"""Database Search (Comet) Results Page."""
import streamlit as st
from pathlib import Path
from src.common.common import page_setup
from src.common.results_helpers import get_workflow_dir
from openms_insight import Table, Heatmap, LinePlot, SequenceView, StateManager

params = page_setup()
st.title("Database Search Results")

st.markdown(
    """
View peptide-spectrum matches (PSMs) identified by **Comet** database search.
Click on a PSM to view the annotated spectrum and peptide sequence.
"""
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

workflow_dir = get_workflow_dir(st.session_state["workspace"])
comet_dir = workflow_dir / "results" / "comet_results"
cache_dir = workflow_dir / "results" / "insight_cache"

if not comet_dir.exists():
    st.info("No database search results available yet. Please run the workflow first.")
    st.page_link("content/workflow_run.py", label="Go to Run", icon="🚀")
    st.stop()

comet_files = sorted(comet_dir.glob("*.idXML"))

if not comet_files:
    st.warning("No identification output files found.")
    st.stop()

selected_file = st.selectbox(
    "Select identification result file",
    comet_files,
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
st.markdown("**Next step:** View rescoring results")
st.page_link("content/results_rescoring.py", label="Go to Rescoring", icon="📈")
