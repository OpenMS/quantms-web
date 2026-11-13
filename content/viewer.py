"""Viewer: every result of a run in views users pick and make, several runs side by side."""
from pathlib import Path

import streamlit as st

from src.common import result_export
from src.common.common import page_setup
from src.common.guide import explain_empty_run, needs_results
from src.common.results_helpers import get_runs
from src.layout_manager import DatasetKind, ViewerSpec, experiment_state, show_viewer
from src.viewer import analyses
from src.viewer.catalog import (
    CATALOG, DEFAULT_LAYOUT, LABELS, MEMBER_HELP, PRESETS, ViewerContext, default_view, member_label,
    migrate_views, renderers,
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
    """What the components of experiment ``ctx`` draw: its run, and the run picked under Files of."""
    return ViewerContext(workspace=workspace, run_id=ctx.dataset, index=ctx.index, params=params,
                         member=ctx.member)


def available(ctx):
    return analyses.available(workspace, ctx.dataset, member=ctx.member)


def add_files_menu(run_id):
    """Add analysis: process more files with the settings of the run shown."""
    with st.popover("Add analysis", icon=":material/library_add:", width="content"):
        if st.button("Process more files", icon=":material/play_circle:", key="viewer_add_files",
                     width="stretch",
                     help="Process more files with these settings and add them to this experiment"):
            analyses.add_files(run_id)


def top_buttons(ctx) -> None:
    """What can be done with the experiment shown: beside its selector, above the views."""
    history_button(runs, ctx.dataset)
    add_files_menu(ctx.dataset)
    result_export.viewer_button(runs, ctx.dataset)


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

show_viewer(ViewerSpec(
    key=KEY, kinds=[RUNS],
    # A run that processed more files for an earlier run is shown with it, as one experiment
    datasets={r["id"]: run_label(r, len(runs.added_runs(r["id"]))) for r in runs.experiments(FINISHED)},
    # Files of: which run's files the components that show one run's files show
    members=runs.experiment, member_mode="files", member_format=lambda m: member_label(workspace, m),
    member_help=MEMBER_HELP,
    label="choose experiment", later_label="Run", buttons=top_buttons, empty=needs_results,
    # A run picked elsewhere (a run that has just finished, or imported from a ZIP) is shown first
    pick=get_selected_run(runs), on_pick=show_picked,
))
