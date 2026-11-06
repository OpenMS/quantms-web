"""The components of the Viewer, the default layout and their renderers.

Component names are stored in saved layouts: never rename one.
"""
from dataclasses import dataclass, field, replace
from importlib import import_module
from pathlib import Path
from typing import Callable, Dict

from src.viewer.scope import experiment

# name, label, group, module in src/viewer/panels
_COMPONENTS = [
    ("run_summary", "Run summary: key figures and files", "Run", "overview"),
    ("sample_counts", "Per-sample counts", "Run", "samples"),
    ("id_quality", "Identification quality", "Run", "qc"),
    ("protein_explorer", "Identified proteins, peptides and spectra", "Identification", "proteins"),
    ("psms", "Peptide-spectrum matches (after FDR)", "Identification", "filtered"),
    ("psms_rescored", "Rescored hits (Percolator)", "Identification", "rescoring"),
    ("psms_search", "Raw search hits", "Identification", "database_search"),
    ("mass_shifts", "Mass shifts", "Identification", "mass_shifts"),
    ("protein_abundance", "Protein and peptide abundances", "Quantification", "abundance"),
    ("group_comparison", "Group comparison: groups, counts and table", "Group comparison", "compare_groups"),
    ("volcano", "Volcano plot", "Group comparison", "volcano"),
    ("pca", "Sample similarity (PCA)", "Group comparison", "pca"),
    ("heatmap", "Heatmap", "Group comparison", "heatmap"),
    ("clustered_heatmap", "Clustered heatmap", "Group comparison", "heatmap_clustered"),
    ("go_enrichment", "GO enrichment", "Group comparison", "enrichment"),
]

# Components that need more than every run has: the analysis each needs (src/viewer/analyses.py)
NEEDS = {
    "run_summary": "identification",
    "sample_counts": "identification",
    "id_quality": "search",
    "protein_explorer": "identification",
    "psms": "identification",
    "psms_rescored": "rescoring",
    "psms_search": "search",
    "mass_shifts": "identification",
    "protein_abundance": "quantification",
    "pca": "quantification",
    "heatmap": ["quantification", "samples"],
    "clustered_heatmap": ["quantification", "samples"],
    "group_comparison": ["quantification", "groups"],
    "volcano": ["quantification", "groups"],
    "go_enrichment": ["quantification", "groups"],
}

# Components computed over the samples of every run of an experiment (a run and
# the runs that processed more files for it). The others show the files of one
# run of it, picked above them.
EXPERIMENT_WIDE = {
    "protein_abundance", "group_comparison", "volcano", "pca", "heatmap",
    "clustered_heatmap", "go_enrichment",
}

CATALOG = [
    {"name": name, "label": label, "group": group, **({"needs": NEEDS[name]} if name in NEEDS else {})}
    for name, label, group, _ in _COMPONENTS
]
MODULES = {name: module for name, _, _, module in _COMPONENTS}
LABELS = {c["name"]: c["label"] for c in CATALOG}

# Ready-made views, one per question users bring to their results, picked with
# the View selector above the Viewer. The first is the default layout.
PRESETS = [
    {
        "name": "overview", "label": "Overview",
        "layout": [["run_summary"], ["id_quality"], ["protein_abundance"]],
        "help": "Did the run work: key figures, identification quality and the protein abundances.",
    },
    {
        "name": "identification", "label": "Identification",
        "layout": [["protein_explorer"], ["psms"], ["psms_rescored"], ["mass_shifts"]],
        "help": "Which proteins and peptides were found, down to the annotated spectra.",
    },
    {
        "name": "quality_control", "label": "Quality control",
        "layout": [["sample_counts"], ["id_quality"], ["psms_search"], ["pca"]],
        "help": "Whether the samples and identifications can be trusted.",
    },
    {
        "name": "differential", "label": "Differential abundance",
        "layout": [["group_comparison"], ["volcano"], ["heatmap"]],
        "help": "Which proteins differ between your groups.",
    },
    {
        "name": "pathways", "label": "Pathways",
        "layout": [["group_comparison"], ["go_enrichment"]],
        "help": "What the differing proteins have in common (GO enrichment).",
    },
]
DEFAULT_LAYOUT = PRESETS[0]["layout"]

# The steps of the group comparison, set in the Edit comparison dialog
COMPARISON_STEPS = [
    ("filtering", "1. Filter proteins"),
    ("imputation", "2. Missing values"),
    ("normalization", "3. Normalize"),
    ("statistical", "4. Statistical test"),
]


@dataclass
class ViewerContext:
    """What a component draws: one run of the workspace, as experiment ``index``."""

    workspace: Path
    run_id: str
    index: int = 0
    params: dict = field(default_factory=dict)


def panel(module: str) -> Callable:
    """The ``render(ctx)`` function of ``src/viewer/panels/<module>.py``."""
    return import_module(f"src.viewer.panels.{module}").render


def draw(module: str, ctx: ViewerContext, name: str = "") -> None:
    """Draw one component (or comparison step) for the context's run.

    A component that shows one run's files, in an experiment of several runs,
    shows the run picked above it (the first run unless picked).
    """
    run_id = ctx.run_id
    if name and name not in EXPERIMENT_WIDE:
        run_id = experiment_member(ctx, name)
    with experiment(run_id, ctx.index):
        panel(module)(replace(ctx, run_id=run_id) if run_id != ctx.run_id else ctx)


def _member_key(run_id: str, index: int) -> str:
    return f"viewer-member-{run_id}-e{index}"


def experiment_member(ctx: ViewerContext, name: str = "") -> str:
    """The run of ``ctx.run_id``'s experiment whose files the components that
    show one run's files show (picked with ``member_picker``, else the first)."""
    import streamlit as st

    from src.common.results_helpers import get_runs

    members = get_runs(ctx.workspace).experiment(ctx.run_id)
    picked = st.session_state.get(_member_key(ctx.run_id, ctx.index))
    return picked if picked in members else ctx.run_id


def member_picker(workspace, run_id: str, index: int) -> None:
    """**Files of**: which run of an experiment of several runs the components
    that show one run's files (PSMs, identified proteins, QC) show. Drawn once
    per experiment, next to its run selector; nothing for a run alone."""
    import streamlit as st

    from src.common.results_helpers import get_runs

    runs = get_runs(workspace)
    members = runs.experiment(run_id)
    if len(members) < 2:
        return

    def label(member: str) -> str:
        files = runs.params(member).get("mzML-files") or []
        count = 1 if isinstance(files, str) else len(files)
        return f"{(runs.get(member) or {}).get('name') or member} ({count} {'file' if count == 1 else 'files'})"

    key = _member_key(run_id, index)
    if st.session_state.get(key) not in members:
        st.session_state[key] = members[0]
    st.segmented_control(
        "Files of", members, key=key, format_func=label,
        help="This experiment's files were processed in several runs (Process more files). "
             "Abundances, the group comparison, PCA and the heatmaps cover all of them; "
             "identifications and quality control show the files of the run picked here.",
    )


def renderers(ctx: ViewerContext) -> Dict[str, Callable[[], None]]:
    """Component name -> function drawing it for the context's run."""
    return {name: (lambda module=module, name=name: draw(module, ctx, name)) for name, module in MODULES.items()}
