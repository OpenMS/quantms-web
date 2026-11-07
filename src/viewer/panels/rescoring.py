"""results_rescoring as a Viewer component."""

import json
import streamlit as st
from pathlib import Path
from src.common.results_helpers import get_run_dir, get_run_params
from src.common.spectrum_view import ensure_psm_views, psm_mzml_files, show_selected_psm_spectrum
from src.common.sequence_display import psm_table, sequence_mode_selector
from openms_insight import Heatmap, SequenceView
from src.viewer.scope import needs_results, scoped, shared_state, PSM_PRIVATE


def render(ctx) -> None:
    run_dir = get_run_dir(ctx.workspace)
    run_params = get_run_params(ctx.workspace)
    perc_dir = run_dir / "results" / "percolator_results"
    cache_dir = run_dir / "results" / "insight_cache"
    if not perc_dir.exists():
        return needs_results()
    prediction_summary = run_dir / "results" / "prediction_features" / "summary.json"
    if prediction_summary.exists():
        with st.container(border=True):
            st.markdown("**Rescored with predicted spectra and retention times.**")
            rows = [
                "| File | Instrument | Collision energy | Matches predicted | Median spectral angle, targets | Median spectral angle, decoys |",
                "|---|---|---|---|---|---|",
            ]
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
        return
    selected_file = st.selectbox(
        "Select rescoring result file", perc_files, format_func=lambda x: x.name, key=scoped("rescoring_file")
    )
    cache_id_prefix = selected_file.stem
    # Check if cache exists
    if not ensure_psm_views(run_dir / "results", selected_file, run_params):
        st.warning("Visualization cache not found. Please re-run the workflow.")
        return
    # Initialize state manager for cross-component linking
    state_manager = shared_state("rescored", PSM_PRIVATE)
    # Load components from cache (no data parameter needed)
    heatmap = Heatmap(cache_id=f"heatmap_{cache_id_prefix}", cache_path=str(cache_dir))
    seq_view = SequenceView(cache_id=f"seqview_{cache_id_prefix}", cache_path=str(cache_dir))
    # Render components
    st.subheader("PSM Overview")
    heatmap(state_manager=state_manager, height=350, key=scoped("rescoring_map"))
    st.subheader("PSM Table")
    mode = sequence_mode_selector("rescored")
    psm_table(cache_dir, f"table_{cache_id_prefix}", mode, state_manager)
    st.subheader("Peptide Sequence")
    seq_view(key=scoped(f"seqview_{cache_id_prefix}"), state_manager=state_manager, height=533)
    st.subheader("MS2 spectrum")
    show_selected_psm_spectrum(
        cache_dir,
        selected_file,
        state_manager,
        seqview_key=scoped(f"seqview_{cache_id_prefix}"),
        workflow_params=run_params,
        mzml_files=psm_mzml_files(run_params, run_dir / "results"),
    )
