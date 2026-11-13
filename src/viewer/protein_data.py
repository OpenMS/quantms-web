"""The identified proteins, peptides and PSMs of a run, and the components drawing them.

The Identification components (protein table, its peptides, their PSMs and the
annotated spectrum of the picked PSM) read the same tables, built once here and
cached, and share the experiment's Insight state, so a protein picked in the
table narrows the peptides, a peptide the PSMs, and a PSM picks the spectrum.
"""
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
import streamlit as st

from src.common.contaminants import is_contaminant
from src.common.fractions import fraction_samples
from src.common.insight_tables import show_insight_table
from src.common.protein_table import build_tables, collect_psms, coverage_html, read_fasta
from src.common.result_filters import add_protein_flags, file_conditions, max_replicates, protein_numeric_filters
from src.common.results_helpers import get_run_dir, get_run_params
from src.viewer.scope import scoped, shared_state

# Columns the linked tables filter on: already shown by the selection above them
HIDDEN = {"visible": False}


def decimals(n: int) -> dict:
    return {"formatter": "money", "formatterParams": {"thousand": False, "precision": n}}


@st.cache_data(show_spinner="Reading the identifications...")
def _load_psms(paths: tuple, fasta: str | None, mtimes: tuple, mzml_names: tuple):
    # mtimes only key the cache, so a re-run of the workflow invalidates it
    psms = collect_psms([Path(p) for p in paths], list(mzml_names))
    return psms, (read_fasta(fasta) if fasta else {})


@st.cache_data(show_spinner="Building protein tables...")
def _load_tables(paths: tuple, fasta: str | None, mtimes: tuple, mzml_names: tuple, conditions: tuple, samples: tuple = ()):
    psms, fasta_entries = _load_psms(paths, fasta, mtimes, mzml_names)
    files = sorted(psms["file"].unique())
    proteins, peptides, psms_long, coverage = build_tables(psms, fasta_entries, files, dict(samples) or None)
    if not proteins.empty:
        proteins = add_protein_flags(proteins, psms_long, dict(conditions), dict(samples) or None)
    return proteins, peptides, psms_long, coverage, fasta_entries


@dataclass
class ProteinData:
    proteins: pd.DataFrame
    peptides: pd.DataFrame
    psms: pd.DataFrame
    coverage: dict
    fasta_entries: dict
    run_dir: Path
    run_params: dict
    cache_dir: Path
    source_files: list
    mzml_names: list
    conditions: dict
    sample_of: dict = field(default_factory=dict)
    fasta_found: bool = True

    @property
    def filter_dir(self) -> Path:
        return self.run_dir / "results" / "psm_filter"


def load(ctx) -> ProteinData | None:
    """The shown run's protein, peptide and PSM tables, None without filtered identifications."""
    from src.WorkflowTest import WorkflowTest, demo_group

    run_dir = get_run_dir(ctx.workspace)
    run_params = get_run_params(ctx.workspace)
    filter_dir = run_dir / "results" / "psm_filter"
    idxml_files = sorted(filter_dir.glob("*.idXML")) if filter_dir.exists() else []
    if not idxml_files:
        return None
    fasta_path = None
    if run_params.get("fasta-file"):
        try:
            fasta_path = Path(WorkflowTest().file_manager.get_files([run_params["fasta-file"]])[0])
        except (ValueError, IndexError):
            fasta_path = None
    fasta_found = bool(fasta_path and fasta_path.is_file())
    mzml_names = [Path(f).name for f in run_params.get("mzML-files", []) or []]
    source_files = [*idxml_files, *([fasta_path] if fasta_found else [])]
    # A run whose samples were split into fractions counts per sample, not per file
    fractions = fraction_samples(filter_dir.parent, run_params) or {}
    sample_of = {f: info["sample"] for f, info in fractions.items()}
    # Conditions from Describe experiment, for "found in n replicates"
    conditions = file_conditions(
        mzml_names,
        lambda name: (fractions.get(name) or {}).get("condition")
        or run_params.get(f"mzML-group-{name}")
        or demo_group(name),
    )
    proteins, peptides, psms_long, coverage, fasta_entries = _load_tables(
        tuple(str(p) for p in idxml_files),
        str(fasta_path) if fasta_found else None,
        tuple(p.stat().st_mtime for p in source_files),
        tuple(mzml_names),
        tuple(sorted(conditions.items())),
        tuple(sorted(sample_of.items())),
    )
    return ProteinData(
        proteins, peptides, psms_long, coverage, fasta_entries, run_dir, run_params,
        run_dir / "results" / "insight_cache", source_files, mzml_names, conditions, sample_of, fasta_found,
    )


def protein_count(data: ProteinData) -> int:
    """Identified target proteins without the contaminants."""
    if data.proteins.empty:
        return 0
    return int((~data.proteins["contaminant"].astype(bool)).sum())


def protein_table(data: ProteinData) -> None:
    """The identified proteins; picking one narrows the components below it."""
    proteins = data.proteins
    if not data.fasta_found:
        st.warning("The FASTA database of this run was not found, so coverage and descriptions cannot be shown.")
    if proteins.empty:
        st.info("No target proteins were identified.")
        return
    show_insight_table(
        proteins.assign(not_contaminant=~proteins["contaminant"]),
        name="proteins",
        cache_dir=data.cache_dir,
        source_files=data.source_files,
        extra_key=f"menu-v1|{'per-sample' if data.sample_of else ''}|{sorted(data.conditions.items())}",
        title="proteins",
        height=420,
        state_manager=shared_state(),
        interactivity={"protein": "accession"},
        column_overrides={
            "coverage (%)": decimals(1),
            "contaminant": {"formatter": "tickCross", "sorter": "boolean", "visible": False},
            "not_contaminant": HIDDEN,
            **(
                {"replicates": {"headerTooltip": "Most files of one condition with a PSM of this protein"}}
                if "replicates" in proteins.columns
                else {}
            ),
        },
        toggle_filters=[{"label": "Hide contaminants", "column": "not_contaminant", "default": True}],
        search=[c for c in ("accession", "gene", "description") if c in proteins.columns],
        search_placeholder="Accession, gene or name",
        numeric_filters=protein_numeric_filters(
            max_replicates(data.conditions, data.sample_of) if "replicates" in proteins.columns else 0
        ),
    )


def protein_peptides(data: ProteinData) -> None:
    """Coverage and peptides of the picked protein; picking a peptide narrows the PSMs."""
    state_manager = shared_state()
    selected = state_manager.get_selection("protein")
    show_insight_table(
        data.peptides,
        name="protein_peptides",
        cache_dir=data.cache_dir,
        source_files=data.source_files,
        extra_key="per-sample" if data.sample_of else "",
        title="peptides",
        height=420,
        state_manager=state_manager,
        filters={"protein": "accession"},
        interactivity={"peptide": "peptide"},
        column_overrides={"accession": HIDDEN},
    )
    # The coverage below the table, so the table starts level with the protein table
    if selected in data.coverage and selected in data.fasta_entries:
        description, sequence = data.fasta_entries[selected]
        mask = data.coverage[selected]
        st.caption(f"{selected} · {description} · {len(sequence)} residues · {100 * mask.mean():.1f} % covered")
        st.html(coverage_html(sequence, mask))
    elif not selected:
        st.caption("Pick a protein to see its coverage and only its peptides.")


def protein_psms(data: ProteinData) -> None:
    """PSMs of the picked protein and peptide; picking one shows its spectrum."""
    show_insight_table(
        data.psms,
        name="protein_psms",
        cache_dir=data.cache_dir,
        source_files=data.source_files,
        title="psms",
        height=860,
        state_manager=shared_state(),
        filters={"protein": "accession", "peptide": "peptide"},
        interactivity={"psm": "psm"},
        extra_key="with-psm-id-v3",  # the cache follows the source files; name the newer psm column and layout
        column_overrides={
            "accession": HIDDEN,
            "peptide": HIDDEN,
            "psm": HIDDEN,
            "RT (s)": decimals(1),
            "m/z": decimals(4),
        },
    )


def psm_spectrum(data: ProteinData) -> None:
    """The annotated spectrum of the picked PSM.

    It comes from the viewer cache the workflow built for the PSM's file (as in
    the PSM explorer), pointed at this one PSM.
    """
    from openms_insight import SequenceView, StateManager

    from src.common.spectrum_view import ensure_psm_views, fragment_tolerance, psm_mzml_files, show_spectrum

    selected_psm = shared_state().get_selection("psm")
    psms = data.psms
    match = psms[psms["psm"] == selected_psm] if selected_psm else psms.iloc[0:0]
    if match.empty:
        st.caption("Pick a PSM to see its annotated spectrum.")
        return
    psm = match.iloc[0]
    run, psm_id = str(selected_psm).rsplit(":", 1)
    mzml_names = data.mzml_names
    has_views = psm["file"] in mzml_names and ensure_psm_views(
        data.run_dir / "results", data.filter_dir / f"{run}.idXML", data.run_params
    )
    if not has_views or not (data.cache_dir / f"seqview_{run}").is_dir():
        st.info("No spectrum viewer data for this run. Re-run the analysis to create it.")
        return
    viewer_state = StateManager(scoped("proteins_page_spectrum"))
    viewer_state.set_selection("identification", int(psm_id))
    viewer_state.set_selection("file", mzml_names.index(psm["file"]))
    viewer_state.set_selection("spectrum", int(psm["scan"]))
    SequenceView(cache_id=f"seqview_{run}", cache_path=str(data.cache_dir))(
        key=scoped("proteins_seqview"), state_manager=viewer_state, height=420,
    )
    tolerance, ppm = fragment_tolerance(data.run_params)
    # The input mzML, or the run's own copy of its identified spectra
    mzml_path = psm_mzml_files(data.run_params, data.run_dir / "results")[mzml_names.index(psm["file"])]
    show_spectrum(
        data.cache_dir, run, viewer_state,
        peptide=psm["peptide"], charge=int(psm["charge"]), file_index=mzml_names.index(psm["file"]),
        scan=int(psm["scan"]), seqview_key=scoped("proteins_seqview"), key=scoped("proteins"),
        tolerance=tolerance, ppm=ppm, mzml_path=mzml_path, height=420,
    )
    st.caption(f"{psm['peptide']} · scan {psm['scan']} · {Path(psm['file']).stem}")


__all__ = [
    "ProteinData", "is_contaminant", "load", "protein_count", "protein_peptides", "protein_psms",
    "protein_table", "psm_spectrum",
]
