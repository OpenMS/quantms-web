"""Viewer: every result of a run in views users pick and make, several runs side by side."""
from pathlib import Path

import streamlit as st

from src.common import result_export
from src.common.common import page_setup
from src.common.guide import explain_empty_run, needs_results
from src.common.results_helpers import get_runs
from src.layout_manager import (
    MAX_EXPERIMENTS, experiment_state, render_experiments, request_analysis, resolve_layout, selected_datasets,
    view_bar,
)
from src.viewer import analyses
from src.viewer.catalog import CATALOG, DEFAULT_LAYOUT, LABELS, PRESETS, ViewerContext, member_picker, renderers
from src.viewer.comparison import edit_comparison_button, open_comparison
from src.workflow.RunManager import FINISHED
from src.workflow.run_widgets import get_selected_run, run_label, set_selected_run

KEY = "viewer_run"

params = page_setup()
st.title("Explore results")
workspace = Path(st.session_state["workspace"])
runs = get_runs(workspace)
finished = runs.list_runs(FINISHED)
# A newest run that ended without results is explained: alone when nothing else
# can be shown, folded away above the runs that can
explain_empty_run(workspace, finished[0]["id"] if finished else None)
if not finished:
    needs_results()

# A run that processed more files for an earlier run is shown with it, as one experiment
datasets = {r["id"]: run_label(r, len(runs.added_runs(r["id"]))) for r in runs.experiments(FINISHED)}
# A run picked elsewhere (a run that has just finished, or imported from a ZIP)
# becomes the first experiment's run
selected = get_selected_run(runs)
if st.session_state.get("_viewer_selected") != selected:
    st.session_state[f"{KEY}_experiment_0"] = selected
    st.session_state["_viewer_selected"] = selected

layout_path = workspace / "layouts" / "viewer.json"
VIEW_KEY = "viewer_view"
# The experiment whose component asked for a missing analysis (its note's button)
MISSING_INDEX = "_viewer_missing_experiment"


def missing_options(missing, view):
    """The dialog's choices for the run of the experiment that asked (the first for a view)."""
    index = st.session_state.pop(MISSING_INDEX, 0)
    run_id = selected_datasets(KEY, index + 1)[index]
    if run_id is None:
        return []
    return analyses.options(workspace, run_id, missing)


def run_of(index: int):
    return selected_datasets(KEY, index + 1)[index]


def add_files():
    run_id = selected_datasets(KEY, 1)[0]
    if run_id is not None:
        analyses.add_files(run_id)


def export_run():
    result_export.open_dialog(runs, selected_datasets(KEY, 1)[0] if datasets else None)


def add_missing(missing, component, run_id, index):
    st.session_state[MISSING_INDEX] = index
    request_analysis(VIEW_KEY, missing, component=component)


view_bar(
    layout_path, CATALOG, PRESETS, DEFAULT_LAYOUT, key=VIEW_KEY,
    available=analyses.available(workspace, selected_datasets(KEY, 1)[0] if datasets else None, 0),
    analyses=analyses.ANALYSES, unavailable_options=missing_options,
    actions=[{"id": "compare", "icon": "⚖️", "help": "Edit comparison: how the groups of the first run "
              "are compared", "on_click": open_comparison},
             {"id": "add-files", "icon": "➕", "label": "Add files", "on_click": add_files,
              "help": "Process more files with the first run's settings and add them to its experiment"},
             {"id": "export", "icon": "📤" if result_export.is_local() else "⬇️",
              "on_click": export_run,
              "help": f"{result_export.verb()} the files of the first run"}],
    max_experiments=MAX_EXPERIMENTS, experiment_label="Run", target="the Viewer",
    # A view of the user's own keeps each run's selection and zoom (📌, ↺)
    states=lambda i: experiment_state("viewer", i, run_of(i)),
    datasets=run_of,
    help_text=(
        "Drag a view to order the bar, untick it to keep it under More. Your own views take up to 3 "
        "components per row; ready-made views can be duplicated. Add runs to compare them side by side."
    ),
)
layout, side_by_side = resolve_layout(layout_path, CATALOG, DEFAULT_LAYOUT)


def experiment_renderers(run_id: str, index: int):
    return renderers(ViewerContext(workspace=workspace, run_id=run_id, index=index, params=params))


def buttons(index: int) -> None:
    run_id = selected_datasets(KEY, index + 1)[index]
    if run_id:
        st.session_state[analyses.INDEX_KEY] = index  # what the run has is asked next, for this experiment
        member_picker(workspace, run_id, index)
        # The first run's comparison opens from ⚖️ in the view bar
        edit_comparison_button(ViewerContext(workspace=workspace, run_id=run_id, index=index, params=params),
                               button=index > 0)


picked = render_experiments(
    layout, side_by_side, datasets, experiment_renderers, key=KEY,
    labels=LABELS, experiment_label="Run", buttons=buttons,
    components=CATALOG, available=lambda run_id: analyses.available(workspace, run_id),
    analyses=analyses.ANALYSES, on_missing=add_missing,
)
if picked and picked[0] and picked[0] != selected:
    set_selected_run(runs, picked[0])
    st.session_state["_viewer_selected"] = picked[0]

