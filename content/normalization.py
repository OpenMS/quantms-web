"""Normalization Page."""

from pathlib import Path
import pandas as pd
import polars as pl
import streamlit as st
from src.common.common import page_setup, save_params, show_fig
from src.common.guide import needs_results, next_steps, page_intro
from src.common.postprocessing_plots import sample_distributions
from src.common.insight_tables import show_insight_table
from src.common.results_helpers import (
    clear_downstream_steps,
    get_abundance_data,
    get_id_column,
    get_sample_group_map,
    postprocessing_param,
)
# Import normalization engine functions from openms_insight
from openms_insight.analysis.normalization import (
    normalize_samples,
    scale_data,
    transform_data,
)

STAT_COLUMNS = ["log2FC", "p-value", "p-adj", "stat"]


def strip_stat_columns(df: pd.DataFrame | None) -> pd.DataFrame | None:
    """Keep preprocessing tables intensity-only before statistical analysis."""
    if df is None:
        return None
    return df.drop(columns=[c for c in STAT_COLUMNS if c in df.columns], errors="ignore")

params = page_setup()
st.title("Fine-tune 3 of 4: Normalize")

page_intro(
    'Make samples comparable by correcting for different amounts of material loaded, and put intensities on a log scale.'
    "\n\n" + 'You do not need this page for a standard analysis: **Compare groups** runs all four steps with the recommended settings. Use it to inspect this step or change it; the plots follow your choice.',
    background="""
- **log2 transformation** turns ratios into differences: a difference of 1 means twice as much.
  It also makes the spread of intensities more even, as statistical tests assume.
- **Median normalization** shifts every sample so its median protein intensity matches the
  others. This corrects for one sample simply containing more material. It assumes most
  proteins do not change, which holds for most experiments.
- **Row scaling** (e.g. auto-scaling) puts every protein on the same scale. It is useful for
  clustering but removes the fold changes, so leave it off for the group comparison.
- **Recommended: log2, median, no scaling.**
""",
    background_title='What normalization does',
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

# Load primary database assets
result = get_abundance_data(st.session_state["workspace"])
if result is None:
    needs_results()

pivot_df, expr_df, group_map = result
pivot_df = strip_stat_columns(pivot_df)
id_col = get_id_column(st.session_state["workspace"], pivot_df)
sample_group_map = get_sample_group_map(st.session_state["workspace"], pivot_df, group_map)

filtered_df = strip_stat_columns(st.session_state.get("filtered_df"))
imputed_df = strip_stat_columns(st.session_state.get("imputed_df"))
normalized_df = strip_stat_columns(st.session_state.get("normalized_df"))
if filtered_df is not None:
    st.session_state["filtered_df"] = filtered_df
if imputed_df is not None:
    st.session_state["imputed_df"] = imputed_df
if normalized_df is not None:
    st.session_state["normalized_df"] = normalized_df

# --- STEP 1: Upstream Pipeline Tracker (Fallback Architecture) ---
if (
    "imputed_df" in st.session_state
    and st.session_state["imputed_df"] is not None
):
    base_df = imputed_df
    st.info("Working on the output of step 2 (Missing values).")
elif (
    "filtered_df" in st.session_state
    and st.session_state["filtered_df"] is not None
):
    base_df = filtered_df
    st.warning("Step 2 (Missing values) has not been applied: working on the output of step 1.")
else:
    base_df = pivot_df
    st.warning("Steps 1 and 2 have not been applied: working on the unprocessed table.")

# 2. Extract actual active sample columns dynamically
sample_cols = [
    c for c in base_df.columns if c not in [id_col, "PeptideSequence", "log2FC", "p-value", "p-adj"]
]

with st.expander(f"Input: {base_df.shape[0]:,} proteins in {len(sample_cols)} samples"):
    show_insight_table(base_df, name="normalization-input")

st.markdown("### Steps so far")
st.caption("Data flows in order: Filtering -> Imputation -> Normalization")

step_rows = [
    {
        "Step": "Filtering",
        "Status": "Done" if filtered_df is not None else "Not run",
        "Rows": filtered_df.shape[0] if filtered_df is not None else "-",
        "Cols": filtered_df.shape[1] if filtered_df is not None else "-",
    },
    {
        "Step": "Imputation",
        "Status": "Done" if imputed_df is not None else "Not run",
        "Rows": imputed_df.shape[0] if imputed_df is not None else "-",
        "Cols": imputed_df.shape[1] if imputed_df is not None else "-",
    },
    {
        "Step": "Normalization",
        "Status": "Done" if normalized_df is not None else "Not run",
        "Rows": normalized_df.shape[0] if normalized_df is not None else "-",
        "Cols": normalized_df.shape[1] if normalized_df is not None else "-",
    },
]
st.dataframe(pd.DataFrame(step_rows), hide_index=True, use_container_width=True)

with st.expander("Show step tables", expanded=False):
    if filtered_df is not None:
        st.markdown("#### Filtering output")
        st.dataframe(filtered_df.head(10), use_container_width=True)
    if imputed_df is not None:
        st.markdown("#### Imputation output")
        st.dataframe(imputed_df.head(10), use_container_width=True)
    if normalized_df is not None:
        st.markdown("#### Normalization output")
        st.dataframe(normalized_df.head(10), use_container_width=True)
    if filtered_df is None and imputed_df is None and normalized_df is None:
        st.info("No preprocessing outputs yet. Start from Filtering.")

st.markdown("---")

# --- SECTION 2: Normalization Parameter Configuration ---
st.subheader("Settings")

# Prepare structural Polars metadata DataFrame required by backend functions
metadata_rows = [{"sample_id": s, "group": sample_group_map[s]} for s in sample_cols if s in sample_group_map]
metadata_pl = pl.DataFrame(
    metadata_rows, schema={"sample_id": pl.String, "group": pl.String}
)

st.caption(
    "Recommended: log2 transformation, median normalization and no row scaling. "
    "Statistical Inference reports the difference of group means as log2FC, so "
    "it expects log2-scaled intensities; median normalization shifts each sample "
    "by a constant and is only meaningful on a log scale; row scaling removes the "
    "fold changes the volcano plot shows."
)

col1, col2, col3 = st.columns(3)

with col1:
    st.markdown("### 🧬 1. Mathematical Transformation")
    transform_options = ["None", "log2", "log10", "square_root", "cube_root"]
    transform_strategy = st.selectbox(
        "Select Transformation",
        options=transform_options,
        index=transform_options.index(
            postprocessing_param(params, "postproc-transform", transform_options)
        ),
        key="postproc-transform",
        help="Compress data dynamic range and stabilize heteroscedastic variance profiles.",
    )

with col2:
    st.markdown("### 🧪 2. Sample Normalization")
    norm_options = ["None", "sum", "median", "pqn", "reference_feature", "quantile"]
    norm_strategy = st.selectbox(
        "Select Normalization",
        options=norm_options,
        index=norm_options.index(
            postprocessing_param(params, "postproc-normalization", norm_options)
        ),
        key="postproc-normalization",
        help="Perform column-wise corrections to account for variable sample loading concentrations.",
    )

    # Conditionally display target input field for reference feature matching
    ref_feature_input = None
    if norm_strategy == "reference_feature":
        ref_feature_input = st.text_input(
            "Reference Protein Name (ID)",
            value="",
            placeholder="e.g., P01234 or GAPDH",
            help=f"Enter the exact unique identifier string matching a key inside the '{id_col}' column.",
        )

with col3:
    st.markdown("### 📊 3. Row Scaling")
    scaling_options = ["None", "mean_centering", "auto_scaling", "pareto_scaling", "range_scaling"]
    scaling_strategy = st.selectbox(
        "Select Scaling Mode",
        options=scaling_options,
        index=scaling_options.index(
            postprocessing_param(params, "postproc-scaling", scaling_options)
        ),
        key="postproc-scaling",
        help="Adjust individual feature weights to make low and high abundance proteins comparable.",
    )

save_params(params)


def run_normalization() -> pl.LazyFrame:
    """Chain the selected transformation, sample normalization and row scaling."""
    processing_lazy = transform_data(
        quantification_data=pl.from_pandas(base_df).lazy(),
        metadata=metadata_pl,
        strategy=transform_strategy,
    )
    processing_lazy = normalize_samples(
        quantification_data=processing_lazy,
        metadata=metadata_pl,
        strategy=norm_strategy,
        id_col=id_col,
        reference_feature=ref_feature_input if norm_strategy == "reference_feature" else None,
    )
    return scale_data(
        quantification_data=processing_lazy,
        metadata=metadata_pl,
        strategy=scaling_strategy,
    )


grouped_samples = metadata_pl["sample_id"].to_list()
with st.expander("📈 Help me choose a method", expanded=True):
    st.caption(
        "Each box is one sample's intensity distribution, before (top) and after "
        "(bottom) the current settings. After a good normalization the medians line up; a sample that stays offset may have a loading or acquisition problem."
    )
    if not grouped_samples:
        st.info("Assign sample groups in Configure to see this plot.")
    else:
        try:
            preview_df = run_normalization().collect().to_pandas()
            show_fig(
                sample_distributions(
                    base_df,
                    preview_df,
                    grouped_samples,
                    sample_group_map,
                    after_is_log=transform_strategy in ("log2", "log10"),
                ),
                "normalization-preview",
            )
        except ValueError as preview_err:
            st.info(f"Preview unavailable: {preview_err}")

# --- SECTION 3: Normalization Pipe Sequential Execution ---
st.markdown("<br>", unsafe_allow_html=True)
if st.button("Apply Normalization Pipelines", type="primary"):

    # Validate reference feature selection if active before hitting polars execution layers
    if norm_strategy == "reference_feature" and not ref_feature_input:
        st.error(
            "❌ Validation Error: Please provide a valid Reference Protein Name to use the 'reference_feature' strategy."
        )
        st.stop()

    try:
        processing_lazy = run_normalization()

        # Finalize and collect pipeline query graph optimizations
        normalized_df = strip_stat_columns(processing_lazy.collect().to_pandas())

        # 💾 Save processing checkpoint inside Session State for Downstream (Statistics Block)
        st.session_state["normalized_df"] = normalized_df
        clear_downstream_steps("normalized_df")
        st.session_state["normalized_transform"] = transform_strategy

        st.success("Successfully executed all selected normalization pipelines!")

    except ValueError as val_err:
        # Gracefully handle validation failures raised from the engine layers (e.g., missing reference protein)
        st.error(f"Engine Configuration Error: {str(val_err)}")
    except Exception as e:
        st.error(f"An unexpected pipeline error occurred: {str(e)}")

# Results are rendered from session state, outside the button block: the Insight
# table triggers a rerun the first time it is drawn, which would reset the button
# and wipe anything shown only while it is pressed.
normalized_df = st.session_state.get("normalized_df")
if normalized_df is not None:
    # Display the finalized transformation matrix view
    st.subheader("Normalized Abundance Table")
    show_insight_table(normalized_df, name="normalization-result")

next_steps(
    ('statistics', 'Next: test every protein for a difference between the groups.'),
    ('compare', 'Back to the overview of the comparison.'),
)
