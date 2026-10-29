"""PCA Results Page."""
import pandas as pd
import polars as pl
import streamlit as st
from src.common.common import page_setup
from src.common.guide import ensure_group_comparison, needs_results, next_steps, page_intro, show_comparison_settings
from src.common.insight_tables import plot_cache, show_insight_table
from src.common.results_helpers import (
    get_abundance_data,
    get_active_table,
    get_id_column,
    get_sample_group_map,
    sample_columns,
)
from openms_insight import PCAPlot

params = page_setup()
st.title("Sample similarity (PCA)")

page_intro(
    'Do replicates of the same group look alike, and do the groups differ? Each dot is one sample, placed so that samples with similar protein profiles sit close together.',
    background="""
- **PCA** (principal component analysis) condenses thousands of protein abundances into a few
  axes that capture most of the differences between samples. **PC1** (x-axis) captures the
  largest share, **PC2** the next largest; the percentages below the plot say how much.
- **What you hope to see**: replicates of one group close together, and the groups apart,
  ideally along PC1.
- **Warning signs**: a sample far away from its own group (possible outlier: check it on
  the Run summary page), or samples grouping by measurement day instead of by condition
  (a batch effect).
""",
    background_title='How to read a PCA plot',
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

# 1. Load abundance data (base wide-format table + sample -> group mapping)
result = get_abundance_data(st.session_state["workspace"])
if result is None:
    needs_results()

pivot_df, expr_df, group_map = result
id_col = get_id_column(st.session_state["workspace"], pivot_df)
sample_group_map = get_sample_group_map(st.session_state["workspace"], pivot_df, group_map)

# --- STEP 1: Use the latest downstream step's output (same table as Statistics) ---
# Plots follow the group comparison, run with the recommended settings if needed
ensure_group_comparison(params, show_errors=False)
base_df, _, _ = get_active_table(pivot_df)

# 2. Extract active sample columns and detect unique biological groups
sample_cols = sample_columns(base_df, sample_group_map)
unique_groups = sorted({sample_group_map[s] for s in sample_cols if s in sample_group_map})

if len(sample_cols) < 2:
    st.info("PCA requires at least 2 samples.")
    st.stop()

if len(unique_groups) < 2:
    st.warning(
        "Only one biological group was detected - points will still be plotted, "
        "but group-based coloring requires 2 or more groups."
    )

with st.expander(f"Data used: {base_df.shape[0]:,} proteins in {len(sample_cols)} samples from {len(unique_groups)} groups"):
    show_insight_table(base_df, name="pca-input")

# --- SECTION 2: PCA Configuration ---


expr_df_wide = base_df.set_index(id_col)[sample_cols]
max_available = expr_df_wide.shape[0]

if max_available <= 20:
    top_n = max_available
    st.caption(f"Using all {top_n} proteins for PCA (dataset too small for variance filtering).")
else:
    top_n = st.slider(
        "Proteins used (most variable first)",
        min_value=20,
        max_value=min(5000, max_available),
        value=min(500, max_available),
        step=10,
        key="pca_top_n",
        help=(
            "PCA is computed only on the N proteins with the highest variance "
            "across samples, to reduce noise from low-variance/uninformative features."
        ),
    )

top_proteins = expr_df_wide.var(axis=1).sort_values(ascending=False).head(top_n).index
expr_df_pca = expr_df_wide.loc[top_proteins].reset_index()

if expr_df_pca.shape[0] < 2:
    st.info("Not enough proteins after variance filtering for PCA.")
    st.stop()

# Prepare structural Polars metadata DataFrame required by PCAPlot
metadata_pl = pl.DataFrame(
    [{"sample_id": s, "group": sample_group_map[s]} for s in sample_cols if s in sample_group_map],
    schema={"sample_id": pl.String, "group": pl.String},
)
pca_lazy = pl.from_pandas(expr_df_pca).lazy()

# 3. Initialize the OpenMS-Insight PCAPlot component (computes PCA internally)
try:
    pca_component = PCAPlot(
        **plot_cache("quantms_pca_plot", expr_df_pca),
        data=pca_lazy,
        metadata=metadata_pl,
        sample_id_field="sample_id",
        group_field="group",
        n_components=5,
        title="Sample PCA",
    )
except ValueError as e:
    st.error(f"PCA computation failed: {e}")
    st.stop()

variance_ratio = pca_component.get_variance_ratio()
pc_columns = pca_component.get_pc_columns()

# 4. Let the user pick which component pair to view (no recomputation needed)
col1, col2 = st.columns(2)
with col1:
    pc_x_label = st.selectbox("X-axis component", pc_columns, index=0, key="pca_pc_x")
with col2:
    default_y_index = 1 if len(pc_columns) > 1 else 0
    pc_y_label = st.selectbox("Y-axis component", pc_columns, index=default_y_index, key="pca_pc_y")

pc_x = int(pc_x_label.replace("PC", ""))
pc_y = int(pc_y_label.replace("PC", ""))

# 5. Render the component
state_manager = st.session_state.get("state")
pca_component(state_manager=state_manager, pc_x=pc_x, pc_y=pc_y, height=600)
st.caption(
    "Each point is a sample, placed by its overall protein profile. Samples of one group should "
    "cluster together and groups should separate; an isolated sample is a candidate outlier worth checking before statistics."
)

st.markdown(
    "**Explained variance:** "
    + ", ".join(f"{col} {ratio * 100:.1f}%" for col, ratio in zip(pc_columns, variance_ratio))
)
st.markdown(f"**Proteins used:** {expr_df_pca.shape[0]} (top {top_n} by variance)")

show_comparison_settings()

next_steps(
    ('clustered', 'The same question as a heatmap, with samples clustered by similarity.'),
    ('volcano', 'Proteins that differ between the groups.'),
    ('overview', "Check an outlying sample's identification numbers."),
)
