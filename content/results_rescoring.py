"""Rescoring (Percolator) Results Page."""
import json
import streamlit as st
from pathlib import Path
from src.common.common import page_setup
from src.common.guide import needs_results, next_steps, page_intro
from src.common.results_helpers import run_picker, get_run_dir, get_run_params
from src.common.spectrum_view import ensure_psm_views, psm_mzml_files, show_selected_psm_spectrum
from src.common.sequence_display import psm_table, sequence_mode_selector
from openms_insight import Heatmap, SequenceView, StateManager

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
- The analysis identifies peptides in three steps: the database search (ProSE) proposes the best
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
run_params = get_run_params(st.session_state["workspace"])
perc_dir = run_dir / "results" / "percolator_results"
cache_dir = run_dir / "results" / "insight_cache"

if not perc_dir.exists():
    needs_results()

prediction_summary = run_dir / "results" / "prediction_features" / "summary.json"
if prediction_summary.exists():
    with st.container(border=True):
        st.markdown(
            "**Rescored with predicted spectra and retention times.** AlphaPeptDeep predicted every "
            "match; Percolator also used how well the measured spectrum matches the predicted one "
            "(spectral angle, 1 = identical) and the retention-time difference. Correct matches "
            "score high, decoys low."
        )
        rows = ["| File | Instrument | Collision energy | Matches predicted | Median spectral angle, targets | Median spectral angle, decoys |",
                "|---|---|---|---|---|---|"]
        for s in json.loads(prediction_summary.read_text(encoding="utf-8")):
            targets, decoys = s.get("median_angle_targets"), s.get("median_angle_decoys")
            rows.append(
                f"| {Path(s['file']).stem} | {s['instrument']} | NCE {s['nce']:g} | "
                f"{s['predicted']:,} of {s['psms']:,} | {targets if targets is None else f'{targets:.2f}'} | "
                f"{decoys if decoys is None else f'{decoys:.2f}'} |"
            )
        st.markdown("\n".join(rows))

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
if not ensure_psm_views(run_dir / "results", selected_file, run_params):
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
    workflow_params=run_params, mzml_files=psm_mzml_files(run_params, run_dir / "results"),
)

next_steps(
    ('filtered', 'The final, filtered matches.'),
    ('qc', 'Overall quality of the identifications.'),
)
