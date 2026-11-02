"""Protein Identification Results Page."""
from pathlib import Path

import streamlit as st
from openms_insight import SequenceView, StateManager

from src.common.common import page_setup
from src.common.guide import needs_results, next_steps, page_intro
from src.common.insight_tables import show_insight_table
from src.common.contaminants import is_contaminant
from src.common.filter_bar import protein_filter_bar
from src.common.protein_table import build_tables, collect_psms, coverage_html, read_fasta
from src.common.fractions import fraction_samples
from src.common.result_filters import (
    add_protein_flags, file_conditions, filter_proteins, filter_psms_by_q, max_replicates,
)
from src.common.results_helpers import run_picker, get_run_dir, get_run_params
from src.common.spectrum_view import fragment_tolerance, show_spectrum
from src.WorkflowTest import WorkflowTest, demo_group

params = page_setup()
run_picker()
wf = WorkflowTest()

st.title("Identified proteins")
page_intro(
    "Every protein found in any sample, with the evidence behind it. **Click a protein** to see "
    "which parts of its sequence were covered by peptides; **click a peptide** to list the spectra "
    "that identified it, and **click one of those** to see the spectrum with its fragment ions "
    "annotated. Use the search fields in the table header to look up your protein of interest.",
    background="""
- A protein is identified through its **peptides**, and each peptide through one or more
  **PSMs** (peptide-spectrum matches: a fragmentation spectrum explained by that peptide).
- **Coverage** is the share of the protein sequence covered by identified peptides. 10 to 40 %
  is normal; trypsin produces some peptides too short or too long to be measured.
- More peptides and PSMs mean more confident identification. Proteins identified by a single
  peptide deserve a closer look.
- All identifications passed a 1 % false discovery rate (unless changed on *Describe experiment*).
""",
    background_title="How identification evidence works",
)

run_dir = get_run_dir(st.session_state["workspace"])
run_params = get_run_params(st.session_state["workspace"])
filter_dir = run_dir / "results" / "psm_filter"
cache_dir = run_dir / "results" / "insight_cache"
idxml_files = sorted(filter_dir.glob("*.idXML")) if filter_dir.exists() else []

if not idxml_files:
    needs_results()

fasta_path = None
if run_params.get("fasta-file"):
    try:
        fasta_path = Path(wf.file_manager.get_files([run_params["fasta-file"]])[0])
    except (ValueError, IndexError):
        fasta_path = None
if fasta_path is None or not fasta_path.is_file():
    st.warning("The FASTA database of this run was not found, so coverage and descriptions cannot be shown.")

mzml_names = [Path(f).name for f in run_params.get("mzML-files", []) or []]
source_files = [*idxml_files, *([fasta_path] if fasta_path and fasta_path.is_file() else [])]


@st.cache_data(show_spinner="Reading the identifications...")
def load_psms(paths: tuple[str, ...], fasta: str | None, mtimes: tuple[float, ...], mzml_names: tuple[str, ...]):
    # mtimes only key the cache, so a re-run of the workflow invalidates it
    psms = collect_psms([Path(p) for p in paths], list(mzml_names))
    fasta_entries = read_fasta(fasta) if fasta else {}
    return psms, fasta_entries


# A run whose samples were split into fractions counts per sample, not per file
fractions = fraction_samples(filter_dir.parent, run_params) or {}
sample_of = {f: info["sample"] for f, info in fractions.items()}


@st.cache_data(show_spinner="Building protein tables...")
def load_tables(paths: tuple[str, ...], fasta: str | None, mtimes: tuple[float, ...], mzml_names: tuple[str, ...],
                q_max: float | None, conditions: tuple, samples: tuple = ()):
    psms, fasta_entries = load_psms(paths, fasta, mtimes, mzml_names)
    psms = filter_psms_by_q(psms, q_max)
    files = sorted(psms["file"].unique())
    proteins, peptides, psms_long, coverage = build_tables(psms, fasta_entries, files, dict(samples) or None)
    if not proteins.empty:
        proteins = add_protein_flags(proteins, psms_long, dict(conditions), dict(samples) or None)
    return proteins, peptides, psms_long, coverage, fasta_entries


# Conditions from Describe experiment, for "found in n replicates"
conditions = file_conditions(
    mzml_names,
    lambda name: (fractions.get(name) or {}).get("condition")
    or run_params.get(f"mzML-group-{name}") or demo_group(name),
)
run_fdr = min(float(run_params.get("psm-fdr-percent", 1.0)), 100.0)
query, q_max, min_peptides, min_replicates = protein_filter_bar(
    "protein-filter", run_fdr if run_fdr < 100.0 else None, max_replicates(conditions, sample_of)
)

proteins, peptides, psms_long, coverage, fasta_entries = load_tables(
    tuple(str(p) for p in idxml_files),
    str(fasta_path) if fasta_path and fasta_path.is_file() else None,
    tuple(p.stat().st_mtime for p in source_files),
    tuple(mzml_names),
    q_max,
    tuple(sorted(conditions.items())),
    tuple(sorted(sample_of.items())),
)

if proteins.empty:
    st.info("No target proteins were identified at this q-value.")
    st.stop()

shown = filter_proteins(proteins, query, min_peptides, min_replicates)
shown_psms = psms_long[psms_long["accession"].isin(shown["accession"])]
clean = ~shown["contaminant"]
c1, c2, c3 = st.columns(3)
c1.metric(
    "Proteins", f"{int(clean.sum()):,}",
    help=f"Matching the filters, without the {int((~clean).sum())} contaminants among them. "
         f"{len(proteins):,} proteins in total at this q-value.",
)
c2.metric("Peptides", f"{shown_psms.loc[~shown_psms['accession'].map(is_contaminant), 'peptide'].nunique():,}")
c3.metric(
    "Target PSMs",
    f"{shown_psms[~shown_psms['accession'].map(is_contaminant)].drop_duplicates(['file', 'scan', 'peptide']).shape[0]:,}",
)
if shown.empty:
    st.info("No protein matches these filters. Loosen them above.")

state_manager = StateManager("proteins_page")
# Columns the linked tables filter on: already shown by the selection above them
HIDDEN = {"visible": False}
# A protein the filters hide is no longer selected, so its peptides go too
if state_manager.get_selection("protein") not in set(shown["accession"]):
    state_manager.clear_selection("protein")


def decimals(n: int) -> dict:
    return {"formatter": "money", "formatterParams": {"thousand": False, "precision": n}}


st.subheader("Proteins")
show_insight_table(
    shown.assign(not_contaminant=~shown["contaminant"]),
    name="proteins", cache_dir=cache_dir, source_files=source_files,
    extra_key=f"filters-v2|{'per-sample' if fractions else ''}|{query}|{q_max}|{min_peptides}|{min_replicates}|{sorted(conditions.items())}",
    title="proteins", height=450, state_manager=state_manager,
    interactivity={"protein": "accession"},
    column_overrides={
        "coverage (%)": decimals(1),
        "contaminant": {"formatter": "tickCross", "sorter": "boolean", "visible": False},
        "not_contaminant": HIDDEN,
        **({"replicates": {"headerTooltip": "Most files of one condition with a PSM of this protein"}}
           if "replicates" in shown.columns else {}),
    },
    toggle_filters=[{"label": "Hide contaminants", "column": "not_contaminant", "default": True}],
)

selected = state_manager.get_selection("protein")
if selected in coverage and selected in fasta_entries:
    description, sequence = fasta_entries[selected]
    mask = coverage[selected]
    st.subheader(f"Coverage of {selected}")
    st.caption(
        f"{description} · {len(sequence)} residues · {mask.sum()} covered ({100 * mask.mean():.1f} %). "
        "Covered residues are highlighted; every occurrence of a peptide counts."
    )
    st.html(coverage_html(sequence, mask))

st.subheader("Peptides of the selected protein")
show_insight_table(
    peptides, name="protein_peptides", cache_dir=cache_dir, source_files=source_files,
    extra_key=f"q|{q_max}",
    title="peptides", height=350, state_manager=state_manager,
    filters={"protein": "accession"}, interactivity={"peptide": "peptide"},
    column_overrides={"accession": HIDDEN},
    extra_key="per-sample" if fractions else "",
)

st.subheader("PSMs of the selected peptide")
st.caption("One row per spectrum that identified the peptide. Click a row to see its annotated spectrum below.")
show_insight_table(
    psms_long, name="protein_psms", cache_dir=cache_dir, source_files=source_files,
    title="psms", height=350, state_manager=state_manager,
    filters={"protein": "accession", "peptide": "peptide"},
    interactivity={"psm": "psm"}, layout="fitColumns",
    extra_key=f"with-psm-id-v2|q|{q_max}",  # the cache follows the source files; name the newer psm column and layout
    column_overrides={
        "accession": HIDDEN, "peptide": HIDDEN, "psm": HIDDEN,
        "RT (s)": decimals(1), "m/z": decimals(4),
    },
)

# The annotated spectrum comes from the viewer cache the workflow built for the
# PSM's run (as on the Peptide-spectrum matches page), pointed at this one PSM
selected_psm = state_manager.get_selection("psm")
match = psms_long[psms_long["psm"] == selected_psm] if selected_psm else psms_long.iloc[0:0]
if not match.empty:
    psm = match.iloc[0]
    run, psm_id = str(selected_psm).rsplit(":", 1)
    st.subheader(f"Spectrum of {psm['peptide']} (scan {psm['scan']}, {Path(psm['file']).stem})")
    if not (cache_dir / f"seqview_{run}").is_dir() or psm["file"] not in mzml_names:
        st.info("No spectrum viewer data for this run. Re-run the analysis to create it.")
    else:
        viewer_state = StateManager("proteins_page_spectrum")
        viewer_state.set_selection("identification", int(psm_id))
        viewer_state.set_selection("file", mzml_names.index(psm["file"]))
        viewer_state.set_selection("spectrum", int(psm["scan"]))
        st.caption(
            "Peaks explained by fragment ions of the peptide are coloured: b ions contain its start, "
            "y ions its end. Click a peak to find it in the sequence. Below the measured spectrum, "
            "mirrored, is the spectrum a deep-learning model (AlphaPeptDeep) predicts for this peptide: "
            "a correct match shows the same strong peaks on both sides."
        )
        SequenceView(cache_id=f"seqview_{run}", cache_path=str(cache_dir))(
            key="proteins_seqview", state_manager=viewer_state, height=420,
        )
        tolerance, ppm = fragment_tolerance(run_params)
        mzml_path = next(
            (Path(f) for f in wf.file_manager.get_files(run_params["mzML-files"]) if Path(f).name == psm["file"]), None,
        )
        show_spectrum(
            cache_dir, run, viewer_state, peptide=psm["peptide"], charge=int(psm["charge"]),
            file_index=mzml_names.index(psm["file"]), scan=int(psm["scan"]), seqview_key="proteins_seqview",
            key="proteins", tolerance=tolerance, ppm=ppm, mzml_path=mzml_path, height=420,
        )

next_steps(
    ("compare", "See whether your protein changes between the groups."),
    ("abundance", "Its intensity in every sample."),
)
