"""enrichment as a Viewer component."""

import pandas as pd
import polars as pl
import streamlit as st
from src.common.common import follow_app_theme
from src.common.guide import ensure_group_comparison, show_comparison_settings
from src.common.insight_tables import show_insight_table
from src.common.results_helpers import get_abundance_data, get_id_column
from openms_insight.analysis.enrichment import calculate_go_enrichment
from src.viewer.scope import needs_results, scoped


def render(ctx) -> None:
    params = ctx.params
    if get_abundance_data(ctx.workspace) is None:
        return needs_results()
    final_statistics_report = ensure_group_comparison(params)
    if final_statistics_report is None:
        return
    # --- STEP 2: Preprocessing Mapping Key Configuration ---
    # Identify target identifier columns dynamically
    abundance_result = get_abundance_data(ctx.workspace)
    id_col = get_id_column(ctx.workspace, abundance_result[0]) if abundance_result else "ProteinName"
    if id_col not in final_statistics_report.columns:
        st.error(f"❌ Structural Error: Column '{id_col}' is missing from the active matrix context.")
        return
    # --- SECTION 1: Parameter Setup & Dynamic Cutoff Labels ---
    st.subheader("Which proteins count as changed?")
    # Check if target p-value should be adjusted or raw based on previous selections (Fallback safely to 'p-adj')
    target_p_col = "p-adj" if "p-adj" in final_statistics_report.columns else "p-value"
    p_label = "Maximum adjusted p-value" if target_p_col == "p-adj" else "Maximum p-value"
    ui_go_col1, ui_go_col2 = st.columns(2)
    with ui_go_col1:
        p_cutoff = st.number_input(
            p_label,
            min_value=0.0001,
            max_value=1.0,
            value=0.05,
            step=0.01,
            format="%.4f",
            help="Proteins below this value (and above the fold change cut-off) are tested for enriched terms.",
            key=scoped("go_p_cutoff"),
        )
    with ui_go_col2:
        fc_cutoff = st.number_input(
            "Minimum change (|log2 fold change|)",
            min_value=0.0,
            max_value=10.0,
            value=1.0,
            step=0.1,
            format="%.2f",
            help="1 = at least twice or at most half as abundant. Up- and down-regulated proteins are tested together.",
            key=scoped("go_fc_cutoff"),
        )
    stats_key = int(pd.util.hash_pandas_object(final_statistics_report[["log2FC", "p-value"]], index=False).sum())
    if st.button("Find enriched GO terms", type="primary", key=scoped("run_go_analysis")):
        with st.spinner("Fetching GO annotations from MyGene.info and testing every term..."):
            status, output = calculate_go_enrichment(
                final_report=pl.from_pandas(final_statistics_report),
                id_col=id_col,
                target_p_col=target_p_col,
                p_cutoff=p_cutoff,
                fc_cutoff=fc_cutoff,
            )
        st.session_state["go_enrichment"] = (stats_key, status, output, target_p_col, p_cutoff, fc_cutoff)
    # Results are rendered from session state, outside the button block: the Insight
    # tables trigger a rerun the first time they are drawn.
    go_result = st.session_state.get("go_enrichment")
    # A result from an earlier comparison no longer matches the table above
    if go_result is not None and go_result[0] == stats_key:
        _, status, output, used_p_col, used_p, used_fc = go_result
        if status == "empty_data":
            st.error("The comparison table has no p-values or fold changes to select proteins from.")
        elif status == "insufficient_proteins":
            st.warning(
                f"Only **{output}** proteins pass the cut-offs ({used_p_col} < {used_p:.4f}, "
                f"|log2 fold change| ≥ {used_fc:.2f}), too few to test. Loosen the cut-offs above."
            )
        elif status == "success":
            st.markdown(
                f"Tested **{output['fg_count']}** changed proteins against **{output['bg_count']}** "
                "quantified proteins as background."
            )
            tabs = st.tabs(["Biological process", "Cellular component", "Molecular function"])
            for tab, go_type in zip(tabs, ["BP", "CC", "MF"]):
                with tab:
                    fig = output["categories"][go_type]["fig"]
                    df_go = output["categories"][go_type]["df"]
                    if fig is not None and df_go is not None:
                        st.plotly_chart(follow_app_theme(fig), use_container_width=True)
                        if isinstance(df_go, pl.DataFrame):
                            df_go = df_go.to_pandas()
                        show_insight_table(df_go, name=f"go-{go_type}", height=400)
                    else:
                        st.info("No term of this category is over-represented among the changed proteins.")
    show_comparison_settings()
