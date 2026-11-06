"""results_database_search as a Viewer component."""

import streamlit as st
from src.common.results_helpers import get_run_dir, get_run_params, search_results_dir
from src.common.spectrum_view import ensure_psm_views, psm_mzml_files, show_selected_psm_spectrum
from src.common.sequence_display import psm_table, sequence_mode_selector
from openms_insight import Heatmap, SequenceView, StateManager
from src.viewer.scope import needs_results, scoped


def render(ctx) -> None:
    st.info(
        "**Score: ln(hyperscore)**, how well the matched b and y fragment ions explain the spectrum "
        "(their number and intensity). Higher is better; the table starts with the best matches."
    )
    run_dir = get_run_dir(ctx.workspace)
    run_params = get_run_params(ctx.workspace)
    search_dir = search_results_dir(run_dir / "results")
    cache_dir = run_dir / "results" / "insight_cache"
    if not search_dir.exists():
        return needs_results()
    search_files = sorted(search_dir.glob("*.idXML"))
    if not search_files:
        st.warning("No identification output files found.")
        return
    selected_file = st.selectbox(
        "Select identification result file", search_files, format_func=lambda x: x.name, key=scoped("search_file")
    )
    cache_id_prefix = selected_file.stem
    # Check if cache exists
    if not ensure_psm_views(run_dir / "results", selected_file, run_params):
        st.warning("Visualization cache not found. Please re-run the workflow.")
        return
    # Initialize state manager for cross-component linking
    state_manager = StateManager(scoped("database_search_state"))
    # Load components from cache (no data parameter needed)
    heatmap = Heatmap(cache_id=f"heatmap_{cache_id_prefix}", cache_path=str(cache_dir))
    seq_view = SequenceView(cache_id=f"seqview_{cache_id_prefix}", cache_path=str(cache_dir))
    # Render components
    st.subheader("PSM Overview")
    heatmap(state_manager=state_manager, height=350, key=scoped("search_map"))
    st.subheader("PSM Table")
    mode = sequence_mode_selector()
    psm_table(cache_dir, f"table_{cache_id_prefix}", mode, state_manager)
    st.subheader("Peptide Sequence")
    seq_view(key=scoped(f"seqview_{cache_id_prefix}"), state_manager=state_manager, height=533)
    st.subheader("MS2 spectrum")
    st.caption(
        "Top: the measured spectrum, with peaks explained by the peptide's fragment ions coloured. Bottom, mirrored: "
        "the spectrum a deep-learning model (AlphaPeptDeep) predicts for the peptide. A correct match shows the same "
        "strong peaks on both sides."
    )
    show_selected_psm_spectrum(
        cache_dir,
        selected_file,
        state_manager,
        seqview_key=f"seqview_{cache_id_prefix}",
        workflow_params=run_params,
        mzml_files=psm_mzml_files(run_params, run_dir / "results"),
    )
