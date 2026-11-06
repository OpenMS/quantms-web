"""The components of the Viewer, the default layout and their renderers.

Component names are stored in saved layouts: never rename one.
"""
from dataclasses import dataclass, field
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

CATALOG = [{"name": name, "label": label, "group": group} for name, label, group, _ in _COMPONENTS]
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


def draw(module: str, ctx: ViewerContext) -> None:
    """Draw one component (or comparison step) for the context's run."""
    with experiment(ctx.run_id, ctx.index):
        panel(module)(ctx)


def renderers(ctx: ViewerContext) -> Dict[str, Callable[[], None]]:
    """Component name -> function drawing it for the context's run."""
    return {name: (lambda module=module: draw(module, ctx)) for name, module in MODULES.items()}
