"""Viewer: every result of a run in views users pick and make, several runs side by side."""
from pathlib import Path

import streamlit as st

from src.common.common import page_setup
from src.common.guide import explain_empty_run, needs_results, page_intro
from src.common.results_helpers import get_runs
from src.layout_manager import MAX_EXPERIMENTS, render_experiments, resolve_layout, selected_datasets, view_bar
from src.viewer.catalog import CATALOG, DEFAULT_LAYOUT, LABELS, PRESETS, ViewerContext, renderers
from src.viewer.comparison import edit_comparison_button, open_comparison
from src.viewer.export import export_section
from src.workflow.RunManager import FINISHED
from src.workflow.run_widgets import get_selected_run, run_label, set_selected_run

KEY = "viewer_run"

params = page_setup()
st.title("Viewer")
page_intro(
    "Every result of your analysis in one place. Each view shows the components for one question, "
    "such as quality control or differential abundance; pick the run to show below the views. **+** "
    "makes a view of your own, ✏️ orders the views and adds runs to compare side by side, ⚖️ "
    "changes how the groups are compared. A protein, peptide or spectrum you pick stays picked "
    "when you switch views. Files of every run are under **Export** at the bottom.",
    background="""
- **Run summary**: did the run work? Proteins, peptides, PSMs and the identification rate.
- **Group comparison**, **Volcano plot**, **PCA**, **Heatmap**, **GO enrichment**: which
  proteins differ between your groups, and what they have in common.
- **Identified proteins**: a protein's coverage, its peptides and their annotated spectra.
- **Identification quality**, **Peptide-spectrum matches**, **Mass shifts**: whether the
  identifications can be trusted, down to single spectra.
""",
    background_title="What the components show",
)

workspace = Path(st.session_state["workspace"])
runs = get_runs(workspace)
finished = runs.list_runs(FINISHED)
# A newest run that ended without results is explained: alone when nothing else
# can be shown, folded away above the runs that can
explain_empty_run(workspace, finished[0]["id"] if finished else None)
if not finished:
    needs_results()

datasets = {r["id"]: run_label(r) for r in finished}
# A run picked elsewhere (👁️ under Export, or a run that has just finished)
# becomes the first experiment's run
selected = get_selected_run(runs)
if st.session_state.get("_viewer_selected") != selected:
    st.session_state[f"{KEY}_experiment_0"] = selected
    st.session_state["_viewer_selected"] = selected

layout_path = workspace / "layouts" / "viewer.json"
view_bar(
    layout_path, CATALOG, PRESETS, DEFAULT_LAYOUT, key="viewer_view",
    actions=[{"id": "compare", "icon": "⚖️", "help": "Edit comparison: how the groups of the first run "
              "are compared", "on_click": open_comparison}],
    max_experiments=MAX_EXPERIMENTS, experiment_label="Run", target="the Viewer",
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
        # The first run's comparison opens from ⚖️ in the view bar
        edit_comparison_button(ViewerContext(workspace=workspace, run_id=run_id, index=index, params=params),
                               button=index > 0)


picked = render_experiments(
    layout, side_by_side, datasets, experiment_renderers, key=KEY,
    labels=LABELS, experiment_label="Run", buttons=buttons,
)
if picked and picked[0] and picked[0] != selected:
    set_selected_run(runs, picked[0])
    st.session_state["_viewer_selected"] = picked[0]

st.divider()
export_section(workspace)
