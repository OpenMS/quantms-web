"""Filtered PSMs Results Page."""
import streamlit as st
from pathlib import Path
from src.common.common import page_setup
from src.common.results_helpers import get_workflow_dir
from src.common.spectrum_view import fragment_tolerance, psm_mzml_files, show_selected_psm_spectrum
from src.common.common import show_fig
from src.common.guide import needs_results, next_steps, page_intro
from src.common.fragment_plots import fragment_error_plot, ion_table, ion_table_html, load_psm_fragments
from src.common.sequence_display import psm_table, sequence_mode_selector
from openms_insight import Heatmap, SequenceView, StateManager

params = page_setup()
st.title("Peptide-spectrum matches")

page_intro(
    'The final, filtered peptide identifications of each run: one row per spectrum explained by a peptide. Select a match to check its spectrum, the fragment ions that support it and their mass errors.',
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
    "**Score: q-value**, the lowest false discovery rate at which this match would still be "
    "accepted. 0.01 means it passes a 1 % FDR; lower is more confident."
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

workflow_dir = get_workflow_dir(st.session_state["workspace"])
filter_dir = workflow_dir / "results" / "psm_filter"
cache_dir = workflow_dir / "results" / "insight_cache"

if not filter_dir.exists():
    needs_results()

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
    workflow_dir=workflow_dir, mzml_files=psm_mzml_files(),
)

st.subheader("Fragment ions")
fragments_tsv = workflow_dir / "results" / "exports" / f"{cache_id_prefix}_fragment_ions.tsv"
psm_id = state_manager.get_selection("identification")
if not fragments_tsv.is_file():
    st.info("No fragment ion export for this run. Re-run the workflow to create it.")
elif psm_id is None:
    st.info("Select a PSM in the table to see its fragment mass errors and ion table.")
else:
    fragments = load_psm_fragments(fragments_tsv, int(psm_id))
    if fragments.empty:
        st.info("No fragment ion of this PSM matched a peak within the tolerance.")
    else:
        tolerance, tolerance_ppm = fragment_tolerance(workflow_dir)

        first = fragments.iloc[0]
        st.caption(
            f"{first['sequence']}, charge {int(first['precursor_charge'])}. Matched with the search's fragment "
            f"tolerance (±{tolerance:g} {'ppm' if tolerance_ppm else 'Da'}, shaded)."
        )
        show_fig(fragment_error_plot(fragments, tolerance, tolerance_ppm), f"fragment-errors-{cache_id_prefix}-{psm_id}")
        values, matched = ion_table(first["sequence"], int(first["precursor_charge"]), fragments)
        st.markdown("**b/y ions** (theoretical m/z; matched ions highlighted, neutral losses not counted)")
        st.html(ion_table_html(values, matched))

next_steps(
    ('qc', 'Overall quality of the identifications.'),
    ('proteins', 'The proteins these peptides add up to.'),
)
