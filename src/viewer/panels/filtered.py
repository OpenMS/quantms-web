"""results_filtered as a Viewer component."""

import pandas as pd
import streamlit as st
from src.common.results_helpers import get_run_dir, get_run_params
from src.common.spectrum_view import ensure_psm_views, fragment_tolerance, psm_mzml_files, show_selected_psm_spectrum
from src.common.fragment_plots import ion_table, ion_table_html, load_psm_fragments
from src.viz.identification import fragment_errors
from src.viz.series import show_series
from src.common.precursor_xic_view import show_precursor_xic
from src.common.result_filters import q_suggestions
from src.common.sequence_display import filtered_psm_table, sequence_mode_selector
from src.viewer.all_psms import all_files_hit, all_files_table_id, file_label
from src.viz.identification.psm_table import ALL_FILES_SELECTION
from openms_insight import Heatmap, SequenceView, Table
from src.viewer.scope import needs_results, scoped, shared_state, PSM_PRIVATE, heading, retitle


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
    stem = getattr(ctx, "file", None)
    if stem:
        # Single run: the PSMs of that run's file
        selected_file = next((f for f in filter_files if f.stem in (f"{stem}_filter", f"{stem}_prose_perc_filter")), None)
        if selected_file is None:
            return needs_results("No PSMs of this run after filtering.")
    else:
        if len(filter_files) > 1 and _all_files_chosen(len(filter_files)):
            return _all_files(cache_dir, run_dir, run_params, filter_files)
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
    # Render components
    head = heading("PSM overview")
    retitle(heatmap, head.get("title"))(
        state_manager=state_manager, height=350, key=scoped("filtered_map"), help=head.get("help")
    )
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
    _psm_details(cache_dir, run_dir, run_params, selected_file, state_manager)


ONE_FILE, ALL_FILES = "One file", "All files"


def _all_files_chosen(count: int) -> bool:
    """**Show**: the PSMs of one file (with its overview) or of all files in one table."""
    choice = st.segmented_control(
        "Show", [ONE_FILE, ALL_FILES], default=ONE_FILE, key=scoped("filtered_scope"),
        format_func=lambda c: f"{c} ({count})" if c == ALL_FILES else c,
        help="All files: every file's PSMs in one table, searchable by file, peptide or protein. "
             "One file: a file's PSMs with their retention time and m/z overview.",
    )
    return choice == ALL_FILES


def _all_files(cache_dir, run_dir, run_params, filter_files) -> None:
    """Every file's PSMs in one table; the picked PSM's sequence, spectrum and fragments below it."""
    for idxml in filter_files:
        ensure_psm_views(run_dir / "results", idxml, run_params)
    mode = sequence_mode_selector("filtered")
    with st.spinner("Putting the PSMs of every file together (first time only)..."):
        table_id = all_files_table_id(cache_dir, filter_files, mode)
    if table_id is None:
        st.warning("Visualization cache not found. Please re-run the workflow.")
        return
    state_manager = shared_state("filtered", PSM_PRIVATE)
    head = heading("PSMs of every file")
    retitle(Table(cache_id=table_id, cache_path=str(cache_dir)), head.get("title"))(
        state_manager=state_manager, height=533, key=scoped("filtered_all_table"),
    )
    hit = all_files_hit(cache_dir, table_id, state_manager.get_selection(ALL_FILES_SELECTION))
    by_stem = {f.stem: f for f in filter_files}
    if hit is None or hit["source"] not in by_stem:
        st.info("Select a PSM in the table to see its sequence, spectrum and fragment ions.")
        return
    selected_file = by_stem[hit["source"]]
    # The picked PSM in its own file's identifiers, which the views below and
    # the One file table read; One file then opens on that file
    state_manager.set_selection("identification", int(hit["id_idx"]))
    state_manager.set_selection("file", int(hit["file_index"]))
    state_manager.set_selection("spectrum", int(hit["scan_id"]))
    st.session_state[scoped("filtered_file")] = selected_file
    st.caption(f"From {file_label(selected_file)}.")
    _psm_details(cache_dir, run_dir, run_params, selected_file, state_manager)


def _psm_details(cache_dir, run_dir, run_params, selected_file, state_manager) -> None:
    """The sequence, spectrum, precursor XIC and fragment ions of the PSM picked in ``selected_file``."""
    cache_id_prefix = selected_file.stem
    seq_view = SequenceView(cache_id=f"seqview_{cache_id_prefix}", cache_path=str(cache_dir))
    retitle(seq_view, "Peptide sequence")(key=scoped(f"seqview_{cache_id_prefix}"), state_manager=state_manager, height=533)
    show_selected_psm_spectrum(
        cache_dir,
        selected_file,
        state_manager,
        seqview_key=scoped(f"seqview_{cache_id_prefix}"),
        workflow_params=run_params,
        mzml_files=psm_mzml_files(run_params, run_dir / "results"),
    )
    show_precursor_xic(
        cache_dir,
        selected_file,
        state_manager,
        mzml_files=psm_mzml_files(run_params),
        key=scoped(f"psm_xic_{cache_id_prefix}"),
    )
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
            unit = "ppm" if tolerance_ppm else "Da"
            error = "mass_error_ppm" if tolerance_ppm else "mass_error_da"
            points = pd.DataFrame({
                "observed_mz": fragments["observed_mz"], error: fragments[error], "ion_type": fragments["ion_type"],
                "fragment": [f"{ion}{'+' * int(z)}" for ion, z in zip(fragments["ion"], fragments["charge"])],
            })
            low, high = float(points["observed_mz"].min()), float(points["observed_mz"].max())
            points = pd.concat([points, pd.DataFrame(
                fragment_errors.tolerance_rows(tolerance, low, high, x="observed_mz", y=error))], ignore_index=True)
            show_series(
                points, f"fragment-errors-{cache_id_prefix}", cache_path=cache_dir,
                key=scoped(f"fragment_errors_{cache_id_prefix}"), height=360,
                **fragment_errors.settings(
                    prsm=None, residue=None, x=("observed_mz", "m/z"), y=(error, f"Mass error ({unit})"),
                    title=f"Fragment ions of {first['sequence']} ({int(first['precursor_charge'])}+): "
                    f"{len(fragments)} matched within ±{tolerance:g} {unit}",
                    svg_filename=f"fragment-errors-{cache_id_prefix}-{psm_id}",
                ),
            )
            values, matched = ion_table(first["sequence"], int(first["precursor_charge"]), fragments)
            st.html(ion_table_html(values, matched))
