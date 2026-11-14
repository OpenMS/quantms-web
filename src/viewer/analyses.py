"""What the Viewer offers when a view needs more than the shown run holds.

The heatmaps compare the samples of a run, so they
need two or more; the group comparison (group comparison, volcano plot, GO
enrichment) needs them in at least two groups. A run of one file, or one whose
samples all share a group, cannot fill them: the view bar locks the views
holding them, and their place in a view shows a note. Both open a dialog that
offers one thing: process more files with the settings of this run, a run of
the added files only, which joins this run's experiment, so the comparison is
computed over both.

A run whose results lack what a component draws (no identifications, no raw
search hits, no rescoring, no abundances) locks it the same way. The components
that show one run's files look at the run picked under **Files of**, those
computed over the experiment at every run of it.
"""
from pathlib import Path

import streamlit as st

from src.common.results_helpers import (
    experiment_fractions,
    get_results_dir,
    search_results_dir,
    experiment_runs,
    get_run_params,
    get_runs,
)

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
# What a run's results hold (src/viewer/panels read the same folders)
IDENTIFICATION = "identification"
SEARCH = "search"
RESCORING = "rescoring"
QUANTIFICATION = "quantification"
_RESULTS = {
    IDENTIFICATION: {
        "label": "identified peptides and proteins",
        "help": "The files shown here have no identifications after FDR filtering.",
    },
    SEARCH: {
        "label": "database search hits",
        "help": "The files shown here have no database search hits.",
    },
    RESCORING: {
        "label": "rescored hits",
        "help": "The files shown here were not rescored with Percolator.",
    },
    QUANTIFICATION: {
        "label": "protein abundances",
        "help": "No run of this experiment has protein abundances.",
    },
}
ANALYSES.update(_RESULTS)
# Set by the Viewer before it asks what a run has: the experiment it draws
CONFIGURE_PAGE = Path("content", "workflow_configure.py")


def _experiment_files(workspace, run_id: str) -> list[str]:
    """The input files of every run of ``run_id``'s experiment."""
    return get_runs(workspace).experiment_files(run_id, "mzML-files")


def run_groups(workspace, run_id: str) -> dict:
    """``{sample: group}`` of the samples of a run's experiment with a group,
    from its current settings.

    The rule of ``load_abundance_data``: a fractionated run groups by the
    conditions of its design, label-free runs by the group of each file, TMT
    runs by the group of each channel (``skip`` leaves a channel out).
    """
    params = get_run_params(workspace, run_id)
    if params.get("analysis-mode", "LFQ") == "TMT":
        return {k: v for k, v in params.items()
                if k.startswith("TMT-group-") and v and str(v).lower() != "skip"}
    fractions = experiment_fractions(workspace, run_id, params)
    if fractions:
        return {info["sample"]: info["condition"] for info in fractions.values() if info["condition"]}
    names = [Path(f).name for f in _experiment_files(workspace, run_id)]
    return {n: params[f"mzML-group-{n}"] for n in names if params.get(f"mzML-group-{n}")}


def run_samples(workspace, run_id: str) -> int:
    """The number of samples a run's experiment measured (files, fractions or channels)."""
    params = get_run_params(workspace, run_id)
    if params.get("analysis-mode", "LFQ") == "TMT":
        return max(len([k for k in params if k.startswith("TMT-group-")]), 2)
    fractions = experiment_fractions(workspace, run_id, params)
    if fractions:
        return len({info["sample"] for info in fractions.values()})
    return len(_experiment_files(workspace, run_id))


def _has(directory: Path, pattern: str) -> bool:
    return directory.is_dir() and any(directory.glob(pattern))


def run_results(workspace, run_id: str) -> set:
    """What the results of one run hold (identification, search, rescoring, quantification)."""
    results = get_results_dir(workspace, run_id)
    has = set()
    if _has(results / "psm_filter", "*.idXML"):
        has.add(IDENTIFICATION)
    if _has(search_results_dir(results), "*.idXML"):
        has.add(SEARCH)
    if _has(results / "percolator_results", "*.idXML"):
        has.add(RESCORING)
    if _has(results / "quant_results", "*.csv"):
        has.add(QUANTIFICATION)
    return has


def available(workspace, run_id, member: str | None = None) -> set | None:
    """The analyses run ``run_id`` can show (None: no run).

    Identifications, search hits and rescoring are those of ``member``, the run
    of the experiment picked under **Files of** (the first when None),
    abundances those of any run of the experiment.
    """
    if run_id is None:
        return None
    from src.viewer.catalog import ViewerContext, experiment_member

    member = experiment_member(ViewerContext(workspace=workspace, run_id=run_id, member=member))
    has = run_results(workspace, member) - {QUANTIFICATION}
    if any(QUANTIFICATION in run_results(workspace, r) for r in experiment_runs(workspace, run_id)):
        has.add(QUANTIFICATION)
    if run_samples(workspace, run_id) >= 2:
        has.add(SAMPLES)
    if len(set(run_groups(workspace, run_id).values())) >= 2:
        has.add(GROUPS)
    return has


def _name(runs, run_id: str) -> str:
    return (runs.get(run_id) or {}).get("name") or run_id


def add_files(run_id: str) -> None:
    """Open Describe experiment with the settings of ``run_id``'s experiment,
    to pick files that are processed and added to it."""
    from src.WorkflowTest import WorkflowTest  # the page set up the workspace already
    WorkflowTest().follow_up(run_id, page=CONFIGURE_PAGE, note="to add more files", add_files="mzML-files")


def options(workspace, run_id: str, missing) -> list:
    """The dialog's choice for run ``run_id`` lacking ``missing``: process more files."""
    name = _name(get_runs(workspace), run_id)
    return [{
        "label": "Process more files with these settings",
        "icon": ":material/library_add:",
        "help": (f"Opens Describe experiment with the settings and groups of {name}, to pick the files "
                 "to add and give them groups. Only those are processed; the Viewer then shows them "
                 f"with the files of {name} as one experiment, and compares the groups over all of them."),
        "on_click": lambda: add_files(run_id),
    }]
