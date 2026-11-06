"""What the Viewer offers when a view needs more than the shown run holds.

The heatmaps compare the samples of a run, so they
need two or more; the group comparison (group comparison, volcano plot, GO
enrichment) needs them in at least two groups. A run of one file, or one whose
samples all share a group, cannot fill them: the view bar locks the views
holding them, and their place in a view shows a note. Both open a dialog with
these choices:

- give the samples groups (Describe experiment; groups apply without a new run),
  when the run has two or more samples;
- process more files with the settings of this run (a follow-up run);
- show a run that has two groups, alone or side by side with this one.
"""
from pathlib import Path

import streamlit as st

from src.common.fractions import fraction_samples
from src.common.results_helpers import get_results_dir, get_run_params, get_runs
from src.layout_manager import load_views, save_views, select_view
from src.workflow.RunManager import FINISHED

# What the components in catalog.NEEDS need
GROUPS = "groups"
SAMPLES = "samples"
ANALYSES = {
    SAMPLES: {
        "label": "two or more samples",
        "help": "Heatmaps compare the samples of a run, so they need at least two.",
    },
    GROUPS: {
        "label": "samples in two groups",
        "help": ("Comparing groups needs samples in at least two groups (e.g. control and treated), "
                 "set in **Samples and groups** on the Describe experiment page."),
    },
}
CONFIGURE_PAGE = Path("content", "workflow_configure.py")


def run_groups(workspace, run_id: str) -> dict:
    """``{sample: group}`` of a run's samples with a group, from its current settings.

    The rule of ``load_abundance_data``: a fractionated run groups by the
    conditions of its design, label-free runs by the group of each file, TMT
    runs by the group of each channel (``skip`` leaves a channel out).
    """
    params = get_run_params(workspace, run_id)
    if params.get("analysis-mode", "LFQ") == "TMT":
        return {k: v for k, v in params.items()
                if k.startswith("TMT-group-") and v and str(v).lower() != "skip"}
    fractions = fraction_samples(get_results_dir(workspace, run_id), params)
    if fractions:
        return {info["sample"]: info["condition"] for info in fractions.values() if info["condition"]}
    files = params.get("mzML-files") or []
    files = [files] if isinstance(files, str) else files
    names = [Path(f).name for f in files]
    return {n: params[f"mzML-group-{n}"] for n in names if params.get(f"mzML-group-{n}")}


def run_samples(workspace, run_id: str) -> int:
    """The number of samples a run measured (files, fractions or channels)."""
    params = get_run_params(workspace, run_id)
    if params.get("analysis-mode", "LFQ") == "TMT":
        return max(len([k for k in params if k.startswith("TMT-group-")]), 2)
    fractions = fraction_samples(get_results_dir(workspace, run_id), params)
    if fractions:
        return len({info["sample"] for info in fractions.values()})
    files = params.get("mzML-files") or []
    return 1 if isinstance(files, str) else len(files)


def available(workspace, run_id) -> set | None:
    """The analyses run ``run_id`` can show (None: no run)."""
    if run_id is None:
        return None
    has = set()
    if run_samples(workspace, run_id) >= 2:
        has.add(SAMPLES)
    if len(set(run_groups(workspace, run_id).values())) >= 2:
        has.add(GROUPS)
    return has


def _name(runs, run_id: str) -> str:
    return (runs.get(run_id) or {}).get("name") or run_id


def _show(view_args, selection_key, index, run_id, view, side_by_side) -> None:
    """Show ``run_id`` (with ``view``): in experiment ``index``, or beside it."""
    path, components, presets, default = view_args
    if side_by_side:
        state = load_views(path, components, presets, default)
        state.update(experiments=2, side_by_side=True)
        if view:
            state["current"] = view
        save_views(path, state)
        st.session_state[selection_key(1)] = run_id
        return
    st.session_state[selection_key(index)] = run_id
    if view:
        select_view(path, view, components, presets, default)


def _follow_up(run_id: str) -> None:
    from src.WorkflowTest import WorkflowTest  # the page set up the workspace already
    WorkflowTest().follow_up(run_id, page=CONFIGURE_PAGE, note="to add more files")


def options(workspace, run_id: str, missing, view, *, selection_key, view_args, index: int = 0) -> list:
    """The dialog's choices for run ``run_id`` lacking ``missing`` (samples, groups)."""
    runs = get_runs(workspace)
    name = _name(runs, run_id)
    choices = []
    if GROUPS in missing and run_samples(workspace, run_id) >= 2:
        choices.append({
            "label": "Give the samples groups",
            "icon": ":material/group_work:",
            "help": ("Opens Describe experiment, where **Samples and groups** says which sample belongs "
                     "to which group. Groups apply to this run without running it again."),
            "on_click": lambda: st.switch_page(CONFIGURE_PAGE),
        })
    choices.append({
        "label": "Process more files with these settings",
        "icon": ":material/library_add:",
        "help": (f"Opens Describe experiment with the settings and samples of {name}, to add files and "
                 "give them groups. The new run is kept next to this one."),
        "primary": not choices,
        "on_click": lambda: _follow_up(run_id),
    })
    view_name = view["name"] if view else None
    having = [r["id"] for r in runs.list_runs(FINISHED)
              if r["id"] != run_id and set(missing) <= (available(workspace, r["id"]) or set())]
    for other in having[:1]:
        other_name = _name(runs, other)
        choices.append({
            "label": f"Show {other_name}",
            "icon": ":material/swap_horiz:",
            "help": f"{other_name} has {' and '.join(ANALYSES[m]['label'] for m in missing)}.",
            "primary": False,
            "on_click": lambda other=other: _show(view_args, selection_key, index, other, view_name, False),
        })
        if index == 0:
            choices.append({
                "label": f"Compare side by side with {other_name}",
                "icon": ":material/splitscreen_right:",
                "help": f"Shows {name} and {other_name} next to each other in this view.",
                "primary": False,
                "on_click": lambda other=other: _show(view_args, selection_key, index, other, view_name, True),
            })
    return choices
