"""Viewer: every result of a run in views users pick and make, several runs side by side."""
from pathlib import Path

import streamlit as st

from src.common import result_export
from src.common.common import page_setup
from src.common.guide import explain_empty_run, needs_results
from src.common.results_helpers import get_runs
from src.layout_manager import DatasetKind, ViewerSpec, experiment_state, show_viewer
from src.viewer import analyses, run_files
from src.viewer.catalog import (
    CATALOG, DEFAULT_LAYOUT, FILE_CATALOG, FILE_DEFAULT_LAYOUT, FILE_PRESETS, LABELS, MEMBER_HELP, PRESETS,
    ViewerContext, default_view, migrate_views, renderers,
)
from src.viewer.comparison import edit_comparison_button
from src.viewer.history_ui import history_button
from src.workflow.RunManager import FINISHED
from src.workflow.run_widgets import get_selected_run, run_label, set_selected_run

KEY = "viewer_run"

params = page_setup()
workspace = Path(st.session_state["workspace"])
runs = get_runs(workspace)
finished = runs.list_runs(FINISHED)
# A newest run that ended without results is explained: alone when nothing else
# can be shown, folded away above the runs that can
explain_empty_run(workspace, finished[0]["id"] if finished else None)
if not finished:
    needs_results()

layout_path = workspace / "layouts" / "viewer.json"
# Views saved before the current ready-made views open on the default view
migrate_views(layout_path)


def run_context(ctx) -> ViewerContext:
    """What the components of experiment ``ctx`` draw: its run, alone (Single run) or with the
    runs of its experiment (Comparison)."""
    return ViewerContext(workspace=workspace, run_id=ctx.dataset, index=ctx.index, params=params,
                         member=ctx.member, single=ctx.member is not None)


def run_of(ctx) -> str:
    """The workflow run shown: the experiment's, or the run that processed the file shown."""
    return run_files.split(ctx.dataset)[0] if run_files.is_file(ctx.dataset) else ctx.dataset


def available(ctx):
    return analyses.available(workspace, ctx.dataset, member=ctx.member)


def file_context(ctx) -> ViewerContext:
    """What the components of a Single run draw: one mzML file of a workflow run."""
    run_id, stem = run_files.split(ctx.dataset)
    return ViewerContext(workspace=workspace, run_id=run_id, index=ctx.index, params=params, file=stem)


def file_members(experiment):
    """The LC-MS runs of an experiment: every file of every workflow run in it."""
    return [run_files.file_id(r, s) for r in runs.experiment(experiment) for s in run_files.stems(runs.params(r))]


def file_label(member):
    run_id, stem = run_files.split(member)
    return stem if len(runs.experiment(run_id)) <= 1 else f"{stem} ({(runs.get(run_id) or {}).get('name') or run_id})"


def add_files_menu(run_id):
    """Add analysis: process more files with the settings of the run shown."""
    with st.popover("Add analysis", icon=":material/library_add:", width="content"):
        if st.button("Process more files", icon=":material/play_circle:", key="viewer_add_files",
                     width="stretch",
                     help="Process more files with these settings and add them to this experiment"):
            analyses.add_files(run_id)


def top_buttons(ctx) -> None:
    """What can be done with the experiment shown: beside its selector, above the views."""
    run_id = run_of(ctx)
    history_button(runs, run_id)
    add_files_menu(run_id)
    result_export.viewer_button(runs, run_id)


def show_picked(ctx) -> None:
    """The experiment shown first is the run the other pages show."""
    if ctx.experiment != get_selected_run(runs):
        set_selected_run(runs, ctx.experiment)


RUNS = DatasetKind(
    name="run", components=CATALOG, presets=PRESETS, default_layout=DEFAULT_LAYOUT, views_path=layout_path,
    renderers=lambda ctx: renderers(run_context(ctx)), labels=LABELS,
    available=available, analyses=analyses.ANALYSES,
    unavailable_options=lambda missing, view, ctx: analyses.options(workspace, ctx.dataset, missing),
    # Quantification for samples in groups, else Identification, until a view is picked
    default_view=lambda ctx: default_view(available(ctx)),
    # A view of the user's own keeps each run's selection and zoom (📌, ↺)
    states=lambda ctx: experiment_state("viewer", ctx.index, ctx.dataset),
    # The comparison steps, beside every experiment's selector
    buttons=lambda ctx: edit_comparison_button(run_context(ctx)),
    view_key="viewer_view", experiment_label="Run", target="the Viewer",
    help_text=(
        "Drag a view to order the bar, untick it to keep it under More. Your own views take up to 3 "
        "components per row; ready-made views can be duplicated. Add runs to compare them side by side."
    ),
)

FILES = DatasetKind(
    name="file", components=FILE_CATALOG, presets=FILE_PRESETS, default_layout=FILE_DEFAULT_LAYOUT,
    views_path=workspace / "layouts" / "viewer-file.json",
    renderers=lambda ctx: renderers(file_context(ctx)), labels=LABELS,
    available=lambda ctx: analyses.available(workspace, run_of(ctx)), analyses=analyses.ANALYSES,
    unavailable_options=lambda missing, view, ctx: analyses.options(workspace, run_of(ctx), missing),
    # One run's picked scan and zoom stay when its view changes
    states=lambda ctx: experiment_state("viewer-file", ctx.index, ctx.dataset),
    view_key="viewer_file_view", experiment_label="Run", target="the Viewer",
    help_text="One LC-MS run: its chromatogram, scans, spectra, peak map and PSMs.",
)

show_viewer(ViewerSpec(
    key=KEY, kinds=[RUNS, FILES], kind_of=lambda d: "file" if run_files.is_file(d) else "run",
    # A run that processed more files for an earlier run is shown with it, as one experiment
    datasets={r["id"]: run_label(r, len(runs.added_runs(r["id"]))) for r in runs.experiments(FINISHED)},
    # Single run | Comparison: one LC-MS run (an mzML file) of the experiment explored on its own,
    # or the experiment as a whole (shown first)
    members=file_members, member_mode="switch", member_format=file_label,
    whole=lambda e: bool(file_members(e)), switch_help=MEMBER_HELP, member_select_label="LC-MS run",
    label="choose experiment", later_label="Run", buttons=top_buttons, empty=needs_results,
    # A run picked elsewhere (a run that has just finished, or imported from a ZIP) is shown first
    pick=get_selected_run(runs), on_pick=show_picked,
))
