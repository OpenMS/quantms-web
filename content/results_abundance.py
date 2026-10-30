"""Abundance (ProteomicsLFQ) Results Page."""
import streamlit as st
import pandas as pd
from pathlib import Path
from src.common.common import page_setup
from src.common.guide import needs_results, next_steps, page_intro
from src.common.results_helpers import get_workflow_dir, get_abundance_data, run_picker, get_run_dir, get_run_params
from src.common.insight_tables import show_insight_table
from src.common.contaminants import CONTAMINANT_PREFIX, contaminant_mask

params = page_setup()
run_picker()
st.title("Protein abundances")

page_intro(
    "How much of each protein is in each sample? One row per protein, one column per sample, "
    "with the raw intensities measured by the mass spectrometer.",
    background="""
- **Intensities are relative**, not absolute amounts: they compare one protein across samples,
  not different proteins with each other (some peptides simply fly better in the instrument).
- **Empty cells** mean the protein was not detected in that sample. The group comparison fills
  them in; see *Fine-tune the comparison → 2. Missing values*.
- Values here are **not yet normalized**. For comparisons, use **Compare groups**, which
  log-transforms and normalizes them first.
- **Contaminants** (keratins from skin, trypsin, BSA, ...) are hidden; switch them on to see them.
""",
    background_title="How to read these tables",
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

run_dir = get_run_dir(st.session_state["workspace"])
quant_dir = run_dir / "results" / "quant_results"
workflow_params = get_run_params(st.session_state["workspace"])
analysis_mode = workflow_params.get("analysis-mode", "LFQ")

if not quant_dir.exists():
    needs_results()

csv_files = sorted(quant_dir.glob("*.csv"))

if not csv_files:
    st.info("No CSV files found in the quantification directory.")
    st.stop()

csv_file = csv_files[0]
cache_dir = run_dir / "results" / "insight_cache"
table_sources = [csv_file, get_workflow_dir(st.session_state["workspace"]) / "params.json"]

show_contaminants = st.toggle(
    "Show contaminants",
    value=False,
    key="abundance-show-contaminants",
    help=f"Proteins with {CONTAMINANT_PREFIX} accessions come from the cRAP contaminants added to the "
    "search database. They are hidden here and left out of the downstream analysis.",
)
# The tables differ with the toggle, so it is part of their cache key
contaminant_key = f"contaminants={show_contaminants}"

def render_protein_table(pivot_df, is_lfq=True):
    """Common function to render the protein-level abundance table"""
    pivot_df = pivot_df.copy()
    if is_lfq:
        st.caption("Each protein's intensity is the sum of its peptides' intensities in that sample.")
    else:
        st.caption(
            "One column per label channel, named by its reporter ion. Each value is the protein's "
            "reporter-ion intensity in that channel, from its most intense peptides; n_peptides "
            "counts the peptides behind it."
        )

    if is_lfq:
        # Handle LFQ mode columns (Raw Intensity)
        id_col = "ProteinName"
        exclude_cols = [id_col, "PeptideSequence"]
        sample_cols = [c for c in pivot_df.columns if c not in exclude_cols]

        display_cols = [id_col] + sample_cols + ["PeptideSequence"]
    else:
        # Handle non-LFQ mode columns (Log2-transformed Intensity)
        id_col = "protein"
        id_col = pivot_df.columns[0]
        exclude_cols = [id_col, "n_proteins", "n_peptides", "protein_score"]
        sample_cols = [c for c in pivot_df.columns if c not in exclude_cols and "ratio" not in c.lower()]

        display_cols = [id_col, "n_peptides"] + sample_cols

    # Filter to available columns, then sort and display
    available_cols = [c for c in display_cols if c in pivot_df.columns]
    view_df = pivot_df[available_cols].copy()
    view_df.insert(1, "Contaminant", contaminant_mask(view_df[id_col]).to_numpy())

    show_insight_table(
        view_df,
        name="abundance_protein",
        cache_dir=cache_dir,
        source_files=table_sources,
        # The cache follows the source files, so name the columns shown as well
        extra_key=f"{contaminant_key};{available_cols}",
        title="protein_abundance",
        # Intensities as whole numbers (the table stores float32)
        column_overrides={
            c: {"formatter": "money", "formatterParams": {"thousand": ",", "precision": 0}, "minWidth": 110}
            for c in sample_cols
        },
        **({} if is_lfq else {"layout": "fitColumns"}),
    )

try:
    df = pd.read_csv(csv_file)

    if df.empty:
        st.info("No data found in this file.")
        st.stop()

    result = get_abundance_data(st.session_state["workspace"], include_contaminants=show_contaminants)
    if not show_contaminants:
        all_proteins = get_abundance_data(st.session_state["workspace"], include_contaminants=True)
        if all_proteins is not None and result is not None:
            hidden = len(all_proteins[0]) - len(result[0])
            if hidden:
                st.caption(f"{hidden} contaminant protein(s) hidden.")

    if analysis_mode == "LFQ":
        protein_tab, psm_tab = st.tabs(["Proteins", "Peptides (per spectrum match)"])

        with protein_tab:
            if result is None:
                st.warning("Could not load abundance data. Please run the workflow first.")
                st.stop()
            
            pivot_df, expr_df, group_map = result
            render_protein_table(pivot_df, is_lfq=True)

        with psm_tab:
            st.caption(
                "The peptide measurements behind the protein table: one row per peptide and charge "
                "state, with its protein, sequence and intensity in every sample."
            )
            psm_df = df.copy()
            if "ProteinName" in psm_df.columns:
                is_con = contaminant_mask(psm_df["ProteinName"]).to_numpy()
                psm_df.insert(1, "Contaminant", is_con)
                if not show_contaminants:
                    psm_df = psm_df[~is_con]
            show_insight_table(
                psm_df,
                name="abundance_psm",
                cache_dir=cache_dir,
                source_files=table_sources,
                extra_key=contaminant_key,
                title="psm_quantification",
            )

    else:
        if result is None:
            needs_results()
        render_protein_table(result[0], is_lfq=False)

except Exception as e:
    st.error(f"Failed to load {csv_file.name}: {e}")

next_steps(
    ("compare", "Test which proteins differ between your groups."),
    ("download", "Get this table as a file."),
)
