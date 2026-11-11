"""results_filtered as a Viewer component."""

import streamlit as st
from src.common.results_helpers import get_run_dir, get_run_params
from src.common.spectrum_view import ensure_psm_views, fragment_tolerance, psm_mzml_files, show_selected_psm_spectrum
from src.common.fragment_plots import fragment_error_plot, ion_table, ion_table_html, load_psm_fragments
from src.common.precursor_xic_view import show_precursor_xic
from src.common.psm_curation_ui import verdict_controls, verdict_list
from src.common.result_filters import q_suggestions
from src.common.sequence_display import filtered_psm_table, sequence_mode_selector
from openms_insight import Heatmap, SequenceView
from src.viewer.scope import needs_results, scoped, show_fig, shared_state, PSM_PRIVATE


def render(ctx) -> None:
    run_dir = get_run_dir(ctx.workspace)
    run_params = get_run_params(ctx.workspace)
    filter_dir = run_dir / "results" / "psm_filter"
    cache_dir = run_dir / "results" / "insight_cache"
    if not filter_dir.exists():
        return needs_results()
    filter_files = sorted(filter_dir.glob("*.idXML"))
    if not filter_files:
        st.warning("No filtering output files found.")
        return
    selected_file = st.selectbox(
        "Select filtering result file", filter_files, format_func=lambda x: x.name, key=scoped("filtered_file")
    )
    cache_id_prefix = selected_file.stem
    # Check if cache exists
    if not ensure_psm_views(run_dir / "results", selected_file, run_params):
        st.warning("Visualization cache not found. Please re-run the workflow.")
        return
    # Initialize state manager for cross-component linking
    state_manager = shared_state("filtered", PSM_PRIVATE)
    # Load components from cache (no data parameter needed)
    heatmap = Heatmap(cache_id=f"heatmap_{cache_id_prefix}", cache_path=str(cache_dir))
    seq_view = SequenceView(cache_id=f"seqview_{cache_id_prefix}", cache_path=str(cache_dir))
    # Render components
    st.subheader("PSM Overview")
    heatmap(state_manager=state_manager, height=350, key=scoped("filtered_map"))
    st.subheader("PSM Table")
    run_fdr = min(float(run_params.get("psm-fdr-percent", 1.0)), 100.0)
    mode = sequence_mode_selector("filtered")
    filtered_psm_table(
        cache_dir,
        f"table_{cache_id_prefix}",
        mode,
        state_manager,
        selected_file,
        q_suggestions=q_suggestions(run_fdr if run_fdr < 100.0 else None),
    )
    verdict_list(run_dir / "results", selected_file, state_manager, scoped(f"verdicts_{cache_id_prefix}"),
                 run_name=run_dir.name)
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
    verdict_controls(run_dir / "results", selected_file, state_manager, scoped(f"verdict_{cache_id_prefix}"))
    show_precursor_xic(
        cache_dir,
        selected_file,
        state_manager,
        mzml_files=psm_mzml_files(run_params),
        key=scoped(f"psm_xic_{cache_id_prefix}"),
    )
    st.subheader("Fragment ions")
    fragments_tsv = run_dir / "results" / "exports" / f"{cache_id_prefix}_fragment_ions.tsv"
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
            tolerance, tolerance_ppm = fragment_tolerance(run_params)

            first = fragments.iloc[0]
            st.caption(
                f"{first['sequence']}, charge {int(first['precursor_charge'])}. Matched with the search's fragment "
                f"tolerance (±{tolerance:g} {'ppm' if tolerance_ppm else 'Da'}, shaded)."
            )
            show_fig(
                fragment_error_plot(fragments, tolerance, tolerance_ppm), f"fragment-errors-{cache_id_prefix}-{psm_id}"
            )
            values, matched = ion_table(first["sequence"], int(first["precursor_charge"]), fragments)
            st.markdown("**b/y ions** (theoretical m/z; matched ions highlighted, neutral losses not counted)")
            st.html(ion_table_html(values, matched))
