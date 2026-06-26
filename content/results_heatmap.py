"""Heatmap Results Page."""
import streamlit as st
import numpy as np
import polars as pl
import plotly.express as px
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from src.common.common import page_setup
from src.common.results_helpers import get_abundance_data
from openms_insight import Heatmap

params = page_setup()
st.title("Heatmap")

st.markdown(
    """
Interactive hierarchically clustered heatmap of protein-level abundance (Z-score normalized).
Powered by OpenMS-Insight multi-resolution engine.
"""
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

# 1. Use the refactored get_abundance_data function (returns only pivot_df and group_map)
result = get_abundance_data(st.session_state["workspace"])
if result is None:
    st.info("Abundance data not available. Please run the workflow and configure sample groups first.")
    st.page_link("content/results_abundance.py", label="Go to Abundance", icon="📋")
    st.stop()

pivot_df, group_map = result

if pivot_df.empty:
    st.info("No data available for heatmap.")
    st.stop()

# 2. Compute expr_df directly and derive sample columns internally
# Select only the actual sample columns, excluding metadata fields like ProteinName.
sample_cols = [c for c in pivot_df.columns if c not in ["ProteinName", "PeptideSequence", "log2FC", "p-adj", "stat", "p-value"]]
expr_df = pivot_df.set_index("ProteinName")[sample_cols]

# 3. UI settings (number of top variance proteins)
top_n = st.slider("Number of proteins (Highest Variance)", 20, 200, 50, key="heatmap_top_n")

# 4. Process data (variance selection -> Z-score normalization)
var_series = expr_df.var(axis=1)
top_proteins = var_series.sort_values(ascending=False).head(top_n).index
heatmap_df = expr_df.loc[top_proteins]

# Compute Z-scores and clean missing/invalid values
heatmap_z = heatmap_df.sub(heatmap_df.mean(axis=1), axis=0).div(heatmap_df.std(axis=1), axis=0)
heatmap_z = heatmap_z.replace([np.inf, -np.inf], np.nan).dropna()

if not heatmap_z.empty:
    # 5. Melt and convert data to Polars to satisfy OpenMS-Insight component requirements
    # Restore the ProteinName row index as a column
    heatmap_z_reset = heatmap_z.reset_index()
    
    # Unpivot the wide-format matrix into long-format (X, Y, Intensity)
    melted_df = heatmap_z_reset.melt(
        id_vars=["ProteinName"], 
        value_vars=sample_cols, 
        var_name="Sample", 
        value_name="Z_score"
    )
    
    # Add sample group mapping if available for heatmap categories
    if group_map:
        melted_df["Group"] = melted_df["Sample"].map(group_map)
        
    # Pack the Pandas DataFrame into a Polars LazyFrame
    heatmap_pl_lazy = pl.from_pandas(melted_df).lazy()

    # 6. Initialize the OpenMS-Insight Heatmap component and map attributes
    # Component spec: X axis (Sample), Y axis (ProteinName), color intensity (Z_score)
    heatmap_component = Heatmap(
        cache_id="quantms_protein_heatmap",
        x_column="Sample",
        y_column="ProteinName",
        data=heatmap_pl_lazy,
        intensity_column="Z_score",        # 🔴 이 컬럼 수치로 색상이 칠해져야 합니다.
        title="Protein Abundance Heatmap (Z-score)",
        x_label="Samples",
        y_label="Proteins",
        colorscale="RdBu",                 # Red-Blue 스케일
        reversescale=True,
        log_scale=False,                   # Z-score는 음수가 있으므로 False 유지
        intensity_label="Z-score",         # 범례 제목을 Z-score로 지정
        category_column=None, 
        min_points=10000,                  # 격자가 잘 표현되도록 점 개수 상한을 넉넉히 지정
    )

    # 7. Render the component
    state_manager = st.session_state.get("state")
    heatmap_component(state_manager=state_manager)

else:
    st.warning("Insufficient data to generate the heatmap.")

st.markdown("---")
st.markdown("**Other visualizations:**")
col1, col2 = st.columns(2)
with col1:
    st.page_link("content/results_volcano.py", label="Volcano Plot", icon="🌋")
with col2:
    st.page_link("content/results_pca.py", label="PCA", icon="📊")