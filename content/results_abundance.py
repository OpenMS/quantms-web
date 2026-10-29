"""Abundance (ProteomicsLFQ) Results Page."""
import streamlit as st
import pandas as pd
from pathlib import Path
from src.common.common import page_setup
from src.common.guide import needs_results, next_steps, page_intro
from src.common.results_helpers import get_workflow_dir, get_abundance_data
from src.workflow.ParameterManager import ParameterManager
from src.common.insight_tables import show_insight_table
from src.common.contaminants import CONTAMINANT_PREFIX, contaminant_mask

params = page_setup()
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

workflow_dir = get_workflow_dir(st.session_state["workspace"])
quant_dir = workflow_dir / "results" / "quant_results"
parameter_manager = ParameterManager(workflow_dir, "TOPP Workflow")

workflow_params = parameter_manager.get_parameters_from_json() 
analysis_mode = workflow_params.get("analysis-mode", "LFQ")

if not quant_dir.exists():
    needs_results()

csv_files = sorted(quant_dir.glob("*.csv"))

if not csv_files:
    st.info("No CSV files found in the quantification directory.")
    st.stop()

csv_file = csv_files[0]
cache_dir = workflow_dir / "results" / "insight_cache"
table_sources = [csv_file, workflow_dir / "params.json"]

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
    st.caption("Each protein's intensity is the sum of its peptides' intensities in that sample.")

    if is_lfq:
        # Handle LFQ mode columns (Raw Intensity)
        id_col = "ProteinName"
        exclude_cols = [id_col, "PeptideSequence"]
        sample_cols = [c for c in pivot_df.columns if c not in exclude_cols]

        display_cols = [id_col] + sample_cols + ["PeptideSequence"]
    else:
        # Handle non-LFQ mode columns (Log2-transformed Intensity)
        id_col = "protein"
        exclude_cols = [id_col, "n_proteins", "n_peptides", "protein_score"]
        sample_cols = [c for c in pivot_df.columns if c not in exclude_cols and "ratio" not in c.lower()]

        display_cols = [id_col] + sample_cols

    # Filter to available columns, then sort and display
    available_cols = [c for c in display_cols if c in pivot_df.columns]
    view_df = pivot_df[available_cols].copy()
    view_df.insert(1, "Contaminant", contaminant_mask(view_df[id_col]).to_numpy())

    show_insight_table(
        view_df,
        name="abundance_protein",
        cache_dir=cache_dir,
        source_files=table_sources,
        extra_key=contaminant_key,
        title="protein_abundance",
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
        pre_processing_tab, protein_tab = st.tabs(["Intensity matrix", "Proteins"])

        if result is None:
            st.info("💡 Please run the workflow first to see results.")
            st.stop()

        pivot_df, expr_df, group_map = result

        with pre_processing_tab:
            show_insight_table(
                pivot_df,
                name="abundance_intensity_matrix",
                cache_dir=cache_dir,
                source_files=table_sources,
                extra_key=contaminant_key,
                title="intensity_matrix",
            )

        with protein_tab:
            render_protein_table(pivot_df, is_lfq=False)

except Exception as e:
    st.error(f"Failed to load {csv_file.name}: {e}")

next_steps(
    ("compare", "Test which proteins differ between your groups."),
    ("download", "Get this table as a file."),
)
