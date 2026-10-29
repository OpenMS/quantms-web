"""Filtered PSMs Results Page."""
import streamlit as st
from pathlib import Path
from src.common.common import page_setup
from src.common.results_helpers import get_workflow_dir
from openms_insight import Table, Heatmap, LinePlot, SequenceView, StateManager
from src.common.common import show_fig
from src.common.fragment_plots import fragment_error_plot, ion_table, ion_table_html, load_psm_fragments
from src.workflow.ParameterManager import ParameterManager
from src.common.sequence_display import psm_table, sequence_mode_selector
from openms_insight import Heatmap, LinePlot, SequenceView, StateManager

params = page_setup()
st.title("Filtered PSMs")

st.markdown(
    """
View FDR-controlled peptide identifications after **IDFilter** processing.
These are high-confidence PSMs that passed the specified FDR threshold.
Click on a PSM to view the annotated spectrum and peptide sequence.
"""
)

st.info(
    "**Score:** The q-value represents the minimum false discovery rate (FDR) at which this PSM "
    "would be accepted. Lower values indicate higher confidence identifications."
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

workflow_dir = get_workflow_dir(st.session_state["workspace"])
filter_dir = workflow_dir / "results" / "psm_filter"
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
        # The search's fragment tolerance: user value, else the workflow's default for it
        workflow_params = ParameterManager(workflow_dir, "TOPP Workflow").get_parameters_from_json()
        comet = "CometAdapter-TMT" if workflow_params.get("analysis-mode", "LFQ") == "TMT" else "CometAdapter"
        comet_params = {**workflow_params.get("_defaults", {}).get(comet, {}), **workflow_params.get(comet, {})}
        tolerance = float(comet_params.get("fragment_mass_tolerance", 0.02))
        tolerance_ppm = comet_params.get("fragment_error_units", "Da") != "Da"

        first = fragments.iloc[0]
        st.caption(
            f"{first['sequence']}, charge {int(first['precursor_charge'])}. Matched with the search's fragment "
            f"tolerance (±{tolerance:g} {'ppm' if tolerance_ppm else 'Da'}, shaded)."
        )
        show_fig(fragment_error_plot(fragments, tolerance, tolerance_ppm), f"fragment-errors-{cache_id_prefix}-{psm_id}")
        values, matched = ion_table(first["sequence"], int(first["precursor_charge"]), fragments)
        st.markdown("**b/y ions** (theoretical m/z; matched ions highlighted, neutral losses not counted)")
        st.html(ion_table_html(values, matched))

st.markdown("---")
st.markdown("**Next step:** View abundance quantification")
st.page_link("content/results_abundance.py", label="Go to Abundance", icon="📋")
