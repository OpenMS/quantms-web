"""Rescoring (Percolator) Results Page."""
import streamlit as st
from pathlib import Path
from src.common.common import page_setup
from src.common.guide import needs_results, next_steps, page_intro
from src.common.results_helpers import run_picker, get_run_dir
from src.common.sequence_display import psm_table, sequence_mode_selector
from openms_insight import Heatmap, LinePlot, SequenceView, StateManager

params = page_setup()
run_picker()
st.title("Rescored hits (Percolator)")

page_intro(
    'The search matches after Percolator has re-scored them with machine learning, before filtering. Use it to see how rescoring changed the confidence of a match.',
    background="""
- A **PSM** (peptide-spectrum match) is one fragmentation spectrum explained by one peptide.
- The **plot** shows every PSM by retention time (when the peptide came off the column) and
  mass; **click a PSM** in the table to see its spectrum with the matched fragment ions
  (b ions from the peptide's start, y ions from its end) and its sequence.
- The analysis identifies peptides in three steps: the database search (Comet) proposes the best
  peptide per spectrum, rescoring (Percolator) re-ranks the matches, and filtering keeps those
  within the false discovery rate. This page shows the result of one step.
""",
    background_title='What you are looking at',
)

st.info(
    "**Score: PEP** (posterior error probability), the probability that this match is wrong. "
    "0.01 means a 1 % chance; lower is more confident."
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

run_dir = get_run_dir(st.session_state["workspace"])
perc_dir = run_dir / "results" / "percolator_results"
cache_dir = run_dir / "results" / "insight_cache"

if not perc_dir.exists():
    needs_results()

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
heatmap = Heatmap(cache_id=f"heatmap_{cache_id_prefix}", cache_path=str(cache_dir))
seq_view = SequenceView(cache_id=f"seqview_{cache_id_prefix}", cache_path=str(cache_dir))
line_plot = LinePlot(cache_id=f"lineplot_{cache_id_prefix}", cache_path=str(cache_dir))

# Render components
st.subheader("PSM Overview")
heatmap(state_manager=state_manager, height=350)

st.subheader("PSM Table")
mode = sequence_mode_selector()
psm_table(cache_dir, f"table_{cache_id_prefix}", mode, state_manager)

st.subheader("Peptide Sequence")
seq_view(key=f"seqview_{cache_id_prefix}", state_manager=state_manager, height=533)

st.subheader("MS2 Spectrum")
line_plot(key=f"lineplot_{cache_id_prefix}", state_manager=state_manager, height=450, sequence_view_key=f"seqview_{cache_id_prefix}")

next_steps(
    ('filtered', 'The final, filtered matches.'),
    ('qc', 'Overall quality of the identifications.'),
)
