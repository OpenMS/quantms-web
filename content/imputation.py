"""Imputation Page."""

from pathlib import Path
import pandas as pd
import polars as pl
import streamlit as st
from src.common.common import page_setup, save_params, show_fig
from src.common.guide import needs_results, next_steps, page_intro
from src.common.postprocessing_plots import imputed_value_preview, missingness_vs_intensity
from src.common.insight_tables import show_insight_table
from src.common.results_helpers import (
    clear_downstream_steps,
    get_abundance_data,
    get_id_column,
    get_sample_group_map,
    group_metadata,
    sample_columns,
    postprocessing_param,
)

# Import imputation algorithms from openms_insight engine
from openms_insight.analysis.imputation import impute_mar, impute_smallest_value

STAT_COLUMNS = ["log2FC", "p-value", "p-adj", "stat"]


def strip_stat_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Keep preprocessing tables intensity-only before statistical analysis."""
    return df.drop(columns=[c for c in STAT_COLUMNS if c in df.columns], errors="ignore")

params = page_setup()
st.title("Fine-tune 2 of 4: Missing values")

page_intro(
    'Fill in the values of proteins that were not detected in some samples, so every protein can be tested.'
    "\n\n" + 'You do not need this page for a standard analysis: **Compare groups** runs all four steps with the recommended settings. Use it to inspect this step or change it; the plots follow your choice.',
    background="""
- **Why values are missing.** A protein is usually missing because it was too faint to detect
  in that sample (*missing not at random*, MNAR), not by chance. Label-free data has many such
  gaps; TMT data has fewer, because all channels of a run are measured in the same spectrum.
- **Recommended: MNAR, row.** A missing value becomes the protein's own lowest measured value,
  i.e. "at most as much as the faintest time we saw it". A protein switched off in one group
  then comes out as clearly lower there, as it should.
- *MAR* (missing at random) fills in the group mean or median. It suits data where dropouts are
  random, but cannot fill a protein missing in a whole group.
- The plots below help decide: if proteins with missing values are generally less intense
  (orange curve lower), the recommended MNAR setting is right.
""",
    background_title='What imputation does',
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

# Load base dataset and clean dictionary keys
result = get_abundance_data(st.session_state["workspace"])
if result is None:
    needs_results()

pivot_df, expr_df, group_map = result
pivot_df = strip_stat_columns(pivot_df)
id_col = get_id_column(st.session_state["workspace"], pivot_df)
sample_group_map = get_sample_group_map(st.session_state["workspace"], pivot_df, group_map)

# 1. Pipeline Checkpoint: Fetch upstream filtered data if available, fallback to raw pivot matrix
if "filtered_df" in st.session_state and st.session_state["filtered_df"] is not None:
    base_df = strip_stat_columns(st.session_state["filtered_df"])
    st.session_state["filtered_df"] = base_df
    st.info("Working on the output of step 1 (Filter proteins).")
else:
    base_df = pivot_df
    st.warning("Step 1 (Filter proteins) has not been applied: working on the unfiltered table.")

# 2. Identify actual sample columns dynamically based on the current active matrix
sample_cols = sample_columns(base_df, sample_group_map)

with st.expander(f"Input: {base_df.shape[0]:,} proteins in {len(sample_cols)} samples"):
    show_insight_table(base_df, name="imputation-input")

# --- SECTION 2: Imputation Configuration ---
st.subheader("Settings")

# Build Polars structural metadata DataFrame
metadata_pl = group_metadata(sample_group_map, sample_cols)

# User selection for core missingness assumption strategy
# In DDA data (label-free in particular) most missing values are proteins below the detection
# limit (MNAR), so the default fills them with the protein's own smallest
# observed intensity. The global minimum sits far below most proteins and turns
# a single random dropout into a large false fold change; MAR group imputation
# cannot fill a protein that is missing in a whole group.
impute_options = ["MNAR (Missing Not At Random)", "MAR (Missing At Random)"]
impute_category = st.selectbox(
    "Why are values missing?",
    options=impute_options,
    index=impute_options.index(
        postprocessing_param(params, "postproc-impute-class", impute_options)
    ),
    key="postproc-impute-class",
    help="MAR uses group metrics (Mean/Median). MNAR fills values below the limit of detection. "
    "Recommended: MNAR with row scope.",
)

# Render algorithmic options sub-menus based on the parent selection
strategy_opt = scope_opt = None  # only the selected class's option is shown
if impute_category == "MAR (Missing At Random)":
    st.markdown(
        "**Group Character Imputation**: Fills missing metrics leveraging sample properties belonging to the same group."
    )
    mar_options = ["median", "mean"]
    strategy_opt = st.radio(
        "Mathematical Strategy",
        options=mar_options,
        index=mar_options.index(
            postprocessing_param(params, "postproc-impute-mar-strategy", mar_options)
        ),
        key="postproc-impute-mar-strategy",
        horizontal=True,
    )

elif impute_category == "MNAR (Missing Not At Random)":
    st.markdown(
        "**Smallest Value Imputation**: Replaces missing items with the minimum values detected to reflect technical dropout limits."
    )
    scope_options = ["row", "global"]
    scope_opt = st.radio(
        "Detection Minimum Scope",
        options=scope_options,
        index=scope_options.index(
            postprocessing_param(params, "postproc-impute-mnar-scope", scope_options)
        ),
        key="postproc-impute-mnar-scope",
        horizontal=True,
        help="'row' targets current protein minimum; 'global' searches the entire mass spectrometry matrix profile.",
    )

save_params(params)


def run_imputation() -> pl.LazyFrame:
    """Apply the selected openms_insight imputation to the input table."""
    quant_lazy = pl.from_pandas(base_df).lazy()
    if impute_category == "MAR (Missing At Random)":
        return impute_mar(
            quantification_data=quant_lazy,
            metadata=metadata_pl,
            group_column="group",
            strategy=strategy_opt,
        )
    return impute_smallest_value(
        quantification_data=quant_lazy, metadata=metadata_pl, scope=scope_opt
    )


grouped_samples = metadata_pl["sample_id"].to_list()
with st.expander("📈 Help me choose a method", expanded=True):
    if not grouped_samples:
        st.info("Assign sample groups in Configure to see these plots.")
    else:
        st.caption(
            "Mean intensity of proteins with and without missing values. If the "
            "orange curve sits clearly lower, values fall below the detection limit (choose MNAR); if the curves overlap, dropout is random (MAR fits)."
        )
        show_fig(missingness_vs_intensity(base_df, grouped_samples), "imputation-missingness")
        st.caption(
            "Where the current setting places the filled-in values (red) relative to "
            "observed ones (blue). They should sit at the low edge of the observed values; a separate spike far below them inflates fold changes."
        )
        preview_df = run_imputation().collect().to_pandas()
        show_fig(imputed_value_preview(base_df, preview_df, grouped_samples), "imputation-preview")

# --- SECTION 3: Imputation Execution ---
if st.button("Apply Imputation", type="primary"):
    imputed_lazy = run_imputation()

    # Resolve lazy graph optimization tree and push to display data frame structure
    imputed_df = strip_stat_columns(imputed_lazy.collect().to_pandas())

    # 💾 Save current output into Session State for down-stream processing (Normalization, Statistics)
    st.session_state["imputed_df"] = imputed_df
    clear_downstream_steps("imputed_df")

    st.success(f"Successfully finalized **{impute_category}** imputation step!")

# Results are rendered from session state, outside the button block: the Insight
# table triggers a rerun the first time it is drawn, which would reset the button
# and wipe anything shown only while it is pressed.
imputed_df = st.session_state.get("imputed_df")
if imputed_df is not None:
    st.subheader("Imputed Result Table")
    show_insight_table(imputed_df, name="imputation-result")

next_steps(
    ('normalization', 'Next: make the samples comparable.'),
    ('compare', 'Back to the overview of the comparison.'),
)
