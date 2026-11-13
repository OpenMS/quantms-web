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
    ("dda_summary", "Identification summary", "Summaries", "dda_summary"),
    ("quant_summary", "Quantification summary", "Summaries", "quant_summary"),
    ("files_table", "Files", "Quality control", "files_table"),
    ("sample_counts", "Counts per sample", "Quality control", "samples"),
    ("id_diagnostics", "Files with few identifications", "Quality control", "diagnostics"),
    ("qc_scores", "Score distribution", "Quality control", "qc_scores"),
    ("qc_mass_error", "Precursor mass error", "Quality control", "qc_mass_error"),
    ("protein_table", "Proteins", "Identification", "protein_table"),
    ("protein_peptides", "Peptides of the picked protein", "Identification", "protein_peptides"),
    ("protein_psms", "PSMs of the picked peptide", "Identification", "protein_psms"),
    ("psm_spectrum", "Annotated spectrum of the picked PSM", "Identification", "psm_spectrum"),
    ("psms", "PSM explorer (after FDR)", "Identification", "filtered"),
    ("psms_rescored", "Rescored hits (Percolator)", "Identification", "rescoring"),
    ("psms_search", "Raw search hits", "Identification", "database_search"),
    ("mass_shifts", "Mass shifts", "Identification", "mass_shifts"),
    ("protein_abundance", "Protein and peptide abundances", "Quantification", "abundance"),
    ("abundance_profile", "Abundance of the picked protein per sample", "Quantification", "abundance_profile"),
    ("group_comparison", "Group comparison", "Quantification", "compare_groups"),
    ("volcano", "Volcano plot", "Quantification", "volcano"),
    ("pca", "Sample similarity (PCA)", "Quantification", "pca"),
    ("heatmap", "Heatmap", "Quantification", "heatmap"),
    ("clustered_heatmap", "Clustered heatmap", "Quantification", "heatmap_clustered"),
    ("go_enrichment", "GO enrichment", "Pathways", "enrichment"),
    ("comparison_missing_values", "Missing values", "Quality control", "comparison_missing_values"),
    ("comparison_imputed_values", "Imputed values", "Quality control", "comparison_imputed_values"),
    ("comparison_normalization", "Normalization", "Quality control", "comparison_normalization"),
    ("comparison_pvalues", "p-value distribution", "Quality control", "comparison_pvalues"),
    # Earlier components, kept for views saved with them
    ("run_summary", "Run summary (summary and files)", "Earlier components", "overview"),
    ("id_quality", "Identification quality (scores, mass error, PSMs per file)", "Earlier components", "qc"),
    ("protein_explorer", "Proteins, peptides and spectra in one", "Earlier components", "proteins"),
]

# Components that need more than every run has: the analysis each needs (src/viewer/analyses.py)
NEEDS = {
    "dda_summary": "identification",
    "quant_summary": "quantification",
    "files_table": "identification",
    "run_summary": "identification",
    "sample_counts": "identification",
    "qc_scores": "search",
    "qc_mass_error": "identification",
    "id_quality": "search",
    "id_diagnostics": "identification",
    "protein_table": "identification",
    "protein_peptides": "identification",
    "protein_psms": "identification",
    "psm_spectrum": "identification",
    "protein_explorer": "identification",
    "psms": "identification",
    "psms_rescored": "rescoring",
    "psms_search": "search",
    "mass_shifts": "identification",
    "protein_abundance": "quantification",
    "abundance_profile": "quantification",
    "pca": "quantification",
    "heatmap": ["quantification", "samples"],
    "clustered_heatmap": ["quantification", "samples"],
    "group_comparison": ["quantification", "groups"],
    "volcano": ["quantification", "groups"],
    "go_enrichment": ["quantification", "groups"],
    "comparison_missing_values": ["quantification", "groups"],
    "comparison_imputed_values": ["quantification", "groups"],
    "comparison_normalization": ["quantification", "groups"],
    "comparison_pvalues": ["quantification", "groups"],
}

# Components computed over the samples of every run of an experiment (a run and
# the runs that processed more files for it). The others show the files of one
# run of it, picked above them.
EXPERIMENT_WIDE = {
    "quant_summary", "protein_abundance", "abundance_profile", "group_comparison", "volcano", "pca", "heatmap",
    "clustered_heatmap", "go_enrichment", "comparison_missing_values", "comparison_imputed_values",
    "comparison_normalization", "comparison_pvalues",
}

CATALOG = [
    {"name": name, "label": label, "group": group, **({"needs": NEEDS[name]} if name in NEEDS else {})}
    for name, label, group, _ in _COMPONENTS
]
MODULES = {name: module for name, _, _, module in _COMPONENTS}
LABELS = {c["name"]: c["label"] for c in CATALOG}

# Ready-made views, one per question users bring to their results. Each opens
# with its own summary, then components that narrow what follows them: a
# protein picked in one is picked in all of them. The first is the default.
PRESETS = [
    {
        "name": "quantification", "label": "Quantification",
        "layout": [
            ["quant_summary"],
            ["volcano", "abundance_profile"],
            ["group_comparison"],
            ["pca", "heatmap"],
        ],
        "help": "Which proteins differ between your groups, and how each one varies across the samples.",
    },
    {
        "name": "pathways", "label": "Pathways",
        "layout": [["volcano", "group_comparison"], ["go_enrichment"]],
        "help": "What the differing proteins have in common (GO enrichment).",
    },
    {
        "name": "identification", "label": "Identification",
        "layout": [
            ["dda_summary"],
            ["protein_table", "protein_peptides"],
            # The PSM's sequence view and spectrum need the full width
            ["protein_psms"],
            ["psm_spectrum"],
            ["mass_shifts"],
        ],
        "help": "Which proteins and peptides were found, down to the annotated spectrum of each match.",
    },
    {
        "name": "quality_control", "label": "Quality control",
        "layout": [
            ["files_table"],
            ["sample_counts", "id_diagnostics"],
            ["qc_scores", "qc_mass_error"],
            ["comparison_missing_values", "comparison_normalization"],
            ["comparison_imputed_values", "comparison_pvalues"],
        ],
        "help": "Whether the identifications and the group comparison can be trusted.",
    },
]
DEFAULT_LAYOUT = PRESETS[0]["layout"]
# Views of earlier releases: a saved file naming one opens on the default view
REMOVED_PRESETS = ("overview", "differential", "comparison_checks")


def default_view(available: set | None) -> str:
    """Quantification for a run whose samples are in groups, else Identification."""
    from src.viewer.analyses import GROUPS

    return "quantification" if available and GROUPS in available else "identification"


def migrate_views(path: Path) -> None:
    """Move a views file saved with the earlier ready-made views to the current ones.

    Entries of removed views are dropped (a current one of them is forgotten, so
    the Viewer opens on its default view), and the ready-made views are put in
    the current order, ahead of the user's own. Done once: a file already in
    order is left alone.
    """
    import json

    try:
        saved = json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return
    views = saved.get("views")
    if not isinstance(views, list):
        return
    names = [v.get("name") for v in views if isinstance(v, dict)]
    preset_names = [p["name"] for p in PRESETS]
    stale = [n for n in names if n in REMOVED_PRESETS]
    in_bar = [n for n in names if n in preset_names]
    if not stale and in_bar == [n for n in preset_names if n in in_bar]:
        return
    kept = [v for v in views if isinstance(v, dict) and v.get("name") not in REMOVED_PRESETS]
    presets = sorted((v for v in kept if v.get("name") in preset_names), key=lambda v: preset_names.index(v["name"]))
    saved["views"] = presets + [v for v in kept if v.get("name") not in preset_names]
    if saved.get("current") in REMOVED_PRESETS:
        saved.pop("current", None)
        saved.pop("layout", None)
    Path(path).write_text(json.dumps(saved, indent=2))


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
    if name in LABELS:
        _title(name)
    run_id = ctx.run_id
    if name and name not in EXPERIMENT_WIDE:
        run_id = experiment_member(ctx, name)
    with experiment(run_id, ctx.index):
        panel(module)(replace(ctx, run_id=run_id) if run_id != ctx.run_id else ctx)


def _title(name: str) -> None:
    """The component's title, with a tooltip linking to its section on Understand your Results."""
    import streamlit as st

    from src.viewer.teaching import tooltip

    st.markdown(f"**{LABELS[name]}**", help=tooltip(name))


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
