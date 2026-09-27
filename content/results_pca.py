"""PCA Results Page."""
import pandas as pd
import polars as pl
import streamlit as st
from src.common.common import page_setup
from src.common.results_helpers import (
    get_abundance_data,
    get_active_table,
    get_id_column,
    get_sample_group_map,
    show_pipeline_banner,
    table_digest,
)
from openms_insight import PCAPlot

params = page_setup()
st.title("PCA Analysis")

st.markdown(
    """
Principal Component Analysis (PCA) of protein-level abundance.
Samples are projected onto their principal components and colored by group assignment to visualize clustering.
"""
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

# 1. Load abundance data (base wide-format table + sample -> group mapping)
result = get_abundance_data(st.session_state["workspace"])
if result is None:
    st.info("Abundance data not available. Please run the workflow and configure sample groups first.")
    st.page_link("content/results_abundance.py", label="Go to Abundance", icon="📋")
    st.stop()

pivot_df, expr_df, group_map = result
id_col = get_id_column(st.session_state["workspace"], pivot_df)
sample_group_map = get_sample_group_map(st.session_state["workspace"], pivot_df, group_map)

# --- STEP 1: Use the latest downstream step's output (same table as Statistics) ---
show_pipeline_banner()
base_df, _, _ = get_active_table(pivot_df)

# 2. Extract active sample columns and detect unique biological groups
sample_cols = [
    c for c in base_df.columns
    if c not in [id_col, "PeptideSequence", "log2FC", "p-adj", "stat", "p-value"]
]
unique_groups = sorted({sample_group_map[s] for s in sample_cols if s in sample_group_map})

if len(sample_cols) < 2:
    st.info("PCA requires at least 2 samples.")
    st.stop()

if len(unique_groups) < 2:
    st.warning(
        "Only one biological group was detected - points will still be plotted, "
        "but group-based coloring requires 2 or more groups."
    )

# --- SECTION 1: Active Input Table Preview ---
st.subheader("Input Table Overview")
st.markdown(
    f"Currently analyzing **{base_df.shape[0]}** rows across **{len(sample_cols)}** samples "
    f"belonging to **{len(unique_groups)} groups** ({', '.join(unique_groups)})."
)
st.dataframe(base_df, use_container_width=True)

st.markdown("---")

# --- SECTION 2: PCA Configuration ---
st.subheader("Configure PCA")

expr_df_wide = base_df.set_index(id_col)[sample_cols]
max_available = expr_df_wide.shape[0]

if max_available <= 20:
    top_n = max_available
    st.caption(f"Using all {top_n} proteins for PCA (dataset too small for variance filtering).")
else:
    top_n = st.slider(
        "Number of proteins (Highest Variance)",
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
        cache_id=f"quantms_pca_plot_{table_digest(expr_df_pca)}",
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

st.markdown(
    "**Explained variance:** "
    + ", ".join(f"{col} {ratio * 100:.1f}%" for col, ratio in zip(pc_columns, variance_ratio))
)
st.markdown(f"**Proteins used:** {expr_df_pca.shape[0]} (top {top_n} by variance)")

st.markdown("---")
st.markdown("**Other visualizations:**")
col1, col2 = st.columns(2)
with col1:
    st.page_link("content/results_volcano.py", label="Volcano Plot", icon="🌋")
with col2:
    st.page_link("content/results_heatmap.py", label="Heatmap", icon="🔥")
