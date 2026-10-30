"""Statistical Inference Page."""

from pathlib import Path
import pandas as pd
import polars as pl
import streamlit as st
from src.common.common import page_setup, save_params, show_fig
from src.common.guide import needs_results, next_steps, page_intro
from src.common.postprocessing_plots import pvalue_histogram
from src.common.insight_tables import show_insight_table
from src.common.results_helpers import (
    get_abundance_data,
    get_id_column,
    get_sample_group_map,
    group_metadata,
    sample_columns,
    postprocessing_param,
)
# Import statistics engine functions from openms_insight
from openms_insight.analysis.statistics import calculate_statistical_tests, adjust_fdr_lazy

params = page_setup()
st.title("Fine-tune 4 of 4: Statistical test")

page_intro(
    'Test every protein for a difference between the groups and correct for testing thousands of proteins at once.'
    "\n\n" + 'You do not need this page for a standard analysis: **Compare groups** runs all four steps with the recommended settings. Use it to inspect this step or change it; the plots follow your choice.',
    background="""
- **limma-like (recommended)**: a t-test that borrows information across all proteins to judge
  each protein's variability. It is more reliable than a plain t-test with the few replicates
  typical in proteomics, and works for two or more groups.
- **Welch**: a plain t-test that allows different variability per group (two groups only).
- **Paired**: for matched samples, e.g. the same patient before and after treatment (two groups,
  same number and order of samples).
- **ANOVA**: tests whether any of several groups differs.
- **Multiple-testing correction**: testing 3,000 proteins at p < 0.05 would call about 150
  significant by chance alone. **Benjamini-Hochberg (BH, recommended)** controls the false
  discovery rate: of the proteins with adjusted p ≤ 0.05, about 5 % are expected to be false.
  *Bonferroni* is much stricter and misses real changes.
""",
    background_title='Which test to choose',
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

# Load primary database assets
result = get_abundance_data(st.session_state["workspace"])
if result is None:
    needs_results()

pivot_df, expr_df, group_map = result
id_col = get_id_column(st.session_state["workspace"], pivot_df)
sample_group_map = get_sample_group_map(st.session_state["workspace"], pivot_df, group_map)

# --- STEP 1: Upstream Pipeline Tracker (Fallback Architecture) ---
if (
    "normalized_df" in st.session_state
    and st.session_state["normalized_df"] is not None
):
    base_df = st.session_state["normalized_df"]
    st.info("Working on the output of step 3 (Normalize).")
elif (
    "imputed_df" in st.session_state
    and st.session_state["imputed_df"] is not None
):
    base_df = st.session_state["imputed_df"]
    st.warning("Step 3 (Normalize) has not been applied: working on the output of step 2.")
elif (
    "filtered_df" in st.session_state
    and st.session_state["filtered_df"] is not None
):
    base_df = st.session_state["filtered_df"]
    st.warning("Steps 2 and 3 have not been applied: working on the output of step 1.")
else:
    base_df = pivot_df
    st.warning("Steps 1 to 3 have not been applied: working on the unprocessed table.")

# log2FC below is the difference of group means, which is only a log2 fold
# change when the input is log2-scaled (the Normalization default).
if st.session_state.get("normalized_df") is None:
    st.warning(
        "⚠️ **Input is not log-transformed**: log2FC is computed as the difference of "
        "group means and is only a log2 fold change on log2-scaled data. Run "
        "**Normalization** (default: log2 + median) before testing."
    )
elif st.session_state.get("normalized_transform", "log2") != "log2":
    st.warning(
        f"⚠️ Normalization used the **{st.session_state['normalized_transform']}** "
        "transformation, so the log2FC column is not a log2 fold change. Use log2 "
        "for fold changes that match the volcano plot axis."
    )

# 2. Extract actual active sample columns and detect unique biological groups
sample_cols = sample_columns(base_df, sample_group_map)
metadata_pl = group_metadata(sample_group_map, sample_cols)
unique_groups = sorted(set(metadata_pl["group"].to_list()))
group_count = len(unique_groups)

with st.expander(f"Input: {base_df.shape[0]:,} proteins in {len(sample_cols)} samples from {group_count} groups ({', '.join(unique_groups)})"):
    show_insight_table(base_df, name="statistical-input")

# --- SECTION 2: Dynamic Statistical Parameter Configuration ---
st.subheader("Settings")

col1, col2 = st.columns(2)

with col1:
    st.markdown("**1. Statistical test**")
    
    if group_count < 2:
        st.error("❌ Statistical testing requires at least 2 unique sample groups.")
        st.stop()

    # Route available method options dynamically based on the group count
    if group_count == 2:
        method_options = ["limma_like", "welch", "paired"]
        help_text = "'limma_like' uses Empirical Bayes variance shrinking. 'welch' is for unequal variances. 'paired' is for dependent samples."
    else:
        method_options = ["limma_like", "anova"]
        help_text = "'limma_like' supports multi-group design matrices. 'anova' computes standard row-wise One-way ANOVA."
        
    selected_method = st.selectbox(
        "Test",
        options=method_options,
        index=method_options.index(
            postprocessing_param(params, "postproc-stat-test", method_options)
        ),
        key="postproc-stat-test",
        help=help_text
    )

with col2:
    st.markdown("**2. Correction for testing many proteins**")
    fdr_options = ["BH", "Bonferroni", "None"]
    selected_fdr = st.selectbox(
        "Correction",
        options=fdr_options,
        index=fdr_options.index(
            postprocessing_param(params, "postproc-fdr", fdr_options)
        ),
        key="postproc-fdr",
        help="'BH' (Benjamini-Hochberg) controls False Discovery Rate. 'Bonferroni' is strict Family-Wise Error Rate control."
    )

save_params(params)

# --- SECTION 3: Statistical Query Execution ---
st.markdown("<br>", unsafe_allow_html=True)
if st.button("Run the test", type="primary"):
    
    # Convert active pandas dataframe into polars lazyframe graph
    stats_lazy = pl.from_pandas(base_df).lazy()
    
    try:
        # Execute Chain 1: Calculate core statistics (Adds log2FC, stat, p-value)
        stats_lazy = calculate_statistical_tests(
            quantification_data=stats_lazy,
            metadata=metadata_pl,
            method=selected_method
        )
        
        # Execute Chain 2: Adjust Multiple Testing (Adds p-adj)
        stats_lazy = adjust_fdr_lazy(
            quantification_data=stats_lazy,
            strategy=selected_fdr
        )
        
        # Resolve lazy graph optimization tree and bring back to pandas memory
        statistics_df = stats_lazy.collect().to_pandas()
        
        # 💾 Save processing checkpoint inside Session State for Downstream (e.g., Volcano plot, Volcano/Heatmap UI)
        st.session_state["statistics_df"] = statistics_df
        st.session_state.pop("comparison_settings", None)
        
        st.success(f"Successfully calculated **{selected_method}** test with **{selected_fdr}** FDR correction!")

    except ValueError as val_err:
        st.error(f"The test could not run: {val_err}")
    except Exception as e:
        st.error(f"An unexpected pipeline error occurred: {str(e)}")

# Results are rendered from session state, outside the button block: the Insight
# table triggers a rerun the first time it is drawn, which would reset the button
# and wipe anything shown only while it is pressed.
statistics_df = st.session_state.get("statistics_df")
if statistics_df is not None:
    # Display the finalized statistics table view
    st.subheader("Result")
    st.caption("log2FC: log2 fold change between the groups; p-adj: p-value adjusted for testing many proteins.")
    show_insight_table(statistics_df, name="statistical-result")

# --- SECTION 4: Test diagnostics ---
if st.session_state.get("statistics_df") is not None:
    with st.expander("📈 Check the test result", expanded=True):
        st.caption(
            "Raw p-values of all proteins. A flat histogram with a peak near 0 means "
            "the test fits and some proteins change; a peak near 1 or a U-shape points to non-log input or the wrong test, and no peak at all means few proteins differ."
        )
        show_fig(pvalue_histogram(st.session_state["statistics_df"]), "statistics-pvalue-histogram")

next_steps(
    ('volcano', 'See the result as a volcano plot.'),
    ('compare', 'Back to the overview of the comparison.'),
)
