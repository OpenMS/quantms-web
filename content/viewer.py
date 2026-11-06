"""Viewer: every result of a run, arranged with Edit layout, several runs side by side."""
from pathlib import Path

import streamlit as st

from src.common.common import page_setup
from src.common.guide import explain_empty_run, needs_results, page_intro
from src.common.results_helpers import get_runs
from src.layout_manager import (
    MAX_EXPERIMENTS,
    layout_button,
    layout_presets,
    render_experiments,
    resolve_layout,
    selected_datasets,
)
from src.viewer.catalog import CATALOG, DEFAULT_LAYOUT, LABELS, PRESETS, ViewerContext, renderers
from src.viewer.comparison import edit_comparison_button
from src.viewer.export import export_section
from src.workflow.RunManager import FINISHED
from src.workflow.run_widgets import get_selected_run, run_label, set_selected_run

KEY = "viewer_run"

params = page_setup()
st.title("Viewer")
page_intro(
    "Every result of your analysis in one place. **View** shows the components for one question, "
    "such as quality control or differential abundance; pick the run to show above them. **Edit "
    "layout** arranges your own components and adds runs to compare side by side. **Edit comparison** changes how the groups are compared. Files of every run are under "
    "**Export** at the bottom.",
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
layout_presets(layout_path, CATALOG, PRESETS, DEFAULT_LAYOUT, key="viewer_view")
layout, side_by_side = resolve_layout(layout_path, CATALOG, DEFAULT_LAYOUT)


def experiment_renderers(run_id: str, index: int):
    return renderers(ViewerContext(workspace=workspace, run_id=run_id, index=index, params=params))


def buttons(index: int) -> None:
    run_id = selected_datasets(KEY, index + 1)[index]
    if index == 0:
        layout_button(
            layout_path, CATALOG, DEFAULT_LAYOUT, key="viewer_layout",
            max_experiments=MAX_EXPERIMENTS, experiment_label="Run", target="the Viewer",
            export_filename="opendda-layout.json",
            help_text=(
                "Drag components into the grid (up to 3 per row). Add runs to compare them side by "
                "side or one below the other, then press **Save**."
            ),
        )
    if run_id:
        edit_comparison_button(ViewerContext(workspace=workspace, run_id=run_id, index=index, params=params))


picked = render_experiments(
    layout, side_by_side, datasets, experiment_renderers, key=KEY,
    labels=LABELS, experiment_label="Run", buttons=buttons,
)
if picked and picked[0] and picked[0] != selected:
    set_selected_run(runs, picked[0])
    st.session_state["_viewer_selected"] = picked[0]

st.divider()
export_section(workspace)
