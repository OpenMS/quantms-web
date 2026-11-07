"""results_abundance as a Viewer component."""

import streamlit as st
import pandas as pd
from src.common.results_helpers import get_workflow_dir, get_abundance_data, get_run_dir, get_run_params
from src.common.insight_tables import show_insight_table
from src.common.contaminants import CONTAMINANT_PREFIX, contaminant_mask
from src.common.fractions import fraction_samples, sample_of
from src.viewer.scope import needs_results, scoped


def render(ctx) -> None:
    run_dir = get_run_dir(ctx.workspace)
    quant_dir = run_dir / "results" / "quant_results"
    workflow_params = get_run_params(ctx.workspace)
    analysis_mode = workflow_params.get("analysis-mode", "LFQ")
    if not quant_dir.exists():
        return needs_results()
    csv_files = sorted(quant_dir.glob("*.csv"))
    if not csv_files:
        st.info("No CSV files found in the quantification directory.")
        return
    csv_file = csv_files[0]
    cache_dir = run_dir / "results" / "insight_cache"
    table_sources = [csv_file, get_workflow_dir(ctx.workspace) / "params.json"]
    show_contaminants = st.toggle(
        "Show contaminants",
        value=False,
        key=scoped("abundance-show-contaminants"),
        help=f"Proteins with {CONTAMINANT_PREFIX} accessions come from the cRAP contaminants added to the "
        "search database. They are hidden here and left out of the downstream analysis.",
    )
    # The tables differ with the toggle, so it is part of their cache key
    contaminant_key = f"contaminants={show_contaminants}"

    def render_protein_table(pivot_df, is_lfq=True):
        """Common function to render the protein-level abundance table"""
        pivot_df = pivot_df.copy()
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
            return

        result = get_abundance_data(ctx.workspace, include_contaminants=show_contaminants)
        if not show_contaminants:
            all_proteins = get_abundance_data(ctx.workspace, include_contaminants=True)
            if all_proteins is not None and result is not None:
                hidden = len(all_proteins[0]) - len(result[0])
                if hidden:
                    st.caption(f"{hidden} contaminant protein(s) hidden.")

        if analysis_mode == "LFQ":
            protein_tab, psm_tab = st.tabs(["Proteins", "Peptides (per spectrum match)"])

            with protein_tab:
                if result is None:
                    st.warning("Could not load abundance data. Please run the workflow first.")
                    return

                pivot_df, expr_df, group_map = result
                render_protein_table(pivot_df, is_lfq=True)

            with psm_tab:
                psm_df = df.copy()
                fractions = fraction_samples(run_dir / "results", workflow_params)
                if fractions and "Reference" in psm_df.columns:
                    # Name the sample each fraction file belongs to, as the protein table does
                    psm_df.insert(
                        psm_df.columns.get_loc("Reference") + 1,
                        "Sample",
                        psm_df["Reference"].map(lambda ref: sample_of(ref, fractions)),
                    )
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
                    extra_key=f"{contaminant_key};samples={bool(fractions)}",
                    title="psm_quantification",
                )

        else:
            if result is None:
                return needs_results()
            render_protein_table(result[0], is_lfq=False)

    except Exception as e:
        st.error(f"Failed to load {csv_file.name}: {e}")
