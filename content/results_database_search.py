"""Database Search (Comet) Results Page."""
import streamlit as st
from pathlib import Path
from src.common.common import page_setup
from src.common.guide import needs_results, next_steps, page_intro
from src.common.results_helpers import run_picker, get_run_dir, get_run_params
from src.common.spectrum_view import psm_mzml_files, show_selected_psm_spectrum
from src.common.sequence_display import psm_table, sequence_mode_selector
from openms_insight import Heatmap, SequenceView, StateManager

params = page_setup()
run_picker()
st.title("Raw search hits (Comet)")

page_intro(
    'The best peptide match for every spectrum straight from the database search, before rescoring and filtering. Many of these matches are wrong; this page is for checking the search itself.',
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
    "**Score: e-value**, the number of random matches expected to score this well by chance. "
    "Lower is better; values far below 1 indicate a confident match."
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

run_dir = get_run_dir(st.session_state["workspace"])
run_params = get_run_params(st.session_state["workspace"])
comet_dir = run_dir / "results" / "comet_results"
cache_dir = run_dir / "results" / "insight_cache"

if not comet_dir.exists():
    needs_results()

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
heatmap = Heatmap(cache_id=f"heatmap_{cache_id_prefix}", cache_path=str(cache_dir))
seq_view = SequenceView(cache_id=f"seqview_{cache_id_prefix}", cache_path=str(cache_dir))

# Render components
st.subheader("PSM Overview")
heatmap(state_manager=state_manager, height=350)

st.subheader("PSM Table")
mode = sequence_mode_selector()
psm_table(cache_dir, f"table_{cache_id_prefix}", mode, state_manager)

st.subheader("Peptide Sequence")
seq_view(key=f"seqview_{cache_id_prefix}", state_manager=state_manager, height=533)

st.subheader("MS2 spectrum")
st.caption(
    "Top: the measured spectrum, with peaks explained by the peptide's fragment ions coloured. Bottom, mirrored: "
    "the spectrum a deep-learning model (AlphaPeptDeep) predicts for the peptide. A correct match shows the same "
    "strong peaks on both sides."
)
show_selected_psm_spectrum(
    cache_dir, selected_file, state_manager, seqview_key=f"seqview_{cache_id_prefix}",
    workflow_params=run_params, mzml_files=psm_mzml_files(run_params),
)

next_steps(
    ('rescoring', 'The same matches after rescoring.'),
    ('filtered', 'The final, filtered matches.'),
)
