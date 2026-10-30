"""Clustered Heatmap Results Page."""
import streamlit as st
import numpy as np
import polars as pl
from src.common.common import page_setup
from src.common.insight_tables import plot_cache
from src.common.guide import ensure_group_comparison, needs_results, next_steps, page_intro, show_comparison_settings
from src.common.results_helpers import (
    run_picker,
    get_abundance_data,
    get_active_table,
    get_id_column,
    get_sample_group_map,
    log2_matrix,
)
from openms_insight import ClusteredHeatmap

params = page_setup()
run_picker()
st.title("Clustered heatmap")

page_intro(
    'The most variable proteins (rows) across all samples (columns). Similar samples and similar proteins are placed next to each other, so patterns shared by a group stand out as blocks of color.',
    background="""
- **Colors** show each protein relative to its own average (Z-score): red is above, blue below,
  white average. This makes rare and abundant proteins comparable.
- **Trees (dendrograms)** on the sides join the most similar samples or proteins first. If the
  sample tree splits along the group color bar at the top, the groups differ consistently.
- **Blocks of rows** with the same color pattern are proteins that rise or fall together,
  often members of one complex or pathway.
""",
    background_title='How to read this heatmap',
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

result = get_abundance_data(st.session_state["workspace"])
if result is None:
    needs_results()

pivot_df, expr_df, group_map = result
id_col = get_id_column(st.session_state["workspace"], pivot_df)
sample_group_map = get_sample_group_map(st.session_state["workspace"], pivot_df, group_map)

# Plot the same protein table as PCA and Statistics: the latest downstream
# step's output, on a log2 scale, proteins with any missing value dropped.
# Plots follow the group comparison, run with the recommended settings if needed
ensure_group_comparison(params, show_errors=False)
base_df, _, is_log2 = get_active_table(pivot_df)
# expr_df from the workflow lists exactly the sample columns in both LFQ and TMT mode
sample_cols = [c for c in expr_df.columns if c in base_df.columns]
expr_df = log2_matrix(base_df, id_col, sample_cols, is_log2).dropna()

if expr_df.empty:
    st.info("No data available for heatmap.")
    st.stop()

top_n = st.slider("Number of proteins (Highest Variance)", 10, 200, 30, key="clustered_heatmap_top_n")

var_series = expr_df.var(axis=1)
top_proteins = var_series.sort_values(ascending=False).head(top_n).index
heatmap_df = expr_df.loc[top_proteins]

heatmap_z = heatmap_df.sub(heatmap_df.mean(axis=1), axis=0).div(heatmap_df.std(axis=1), axis=0)
heatmap_z = heatmap_z.replace([np.inf, -np.inf], np.nan).dropna()

if heatmap_z.empty:
    st.warning("Insufficient data to generate the heatmap.")
    st.stop()

heatmap_z_reset = heatmap_z.reset_index()
heatmap_lazy = pl.from_pandas(heatmap_z_reset).lazy()

sample_cols = expr_df.columns.tolist()
metadata_pl = pl.DataFrame(
    [{"sample_id": s, "group": sample_group_map[s]} for s in sample_cols if s in sample_group_map],
    schema={"sample_id": pl.String, "group": pl.String},
)

# Assign group annotation-bar colors in sorted-group order (matching how
# ClusteredHeatmap._preprocess() orders unique groups internally).
group_palette = [
    "#00BFC4",  # teal
    "#F8766D",  # salmon
    "#7CAE00",  # yellow-green
    "#C77CFF",  # lavender purple
    "#E7B800",  # gold/amber
    "#619CFF",  # blue
    "#FF61C3",  # pink/magenta
    "#00BA38",  # green
    "#FF8C42",  # orange
    "#00B0F6",  # sky blue
]
unique_groups = sorted(set(sample_group_map.values()))
group_colors = {g: group_palette[i % len(group_palette)] for i, g in enumerate(unique_groups)}

heatmap_component = ClusteredHeatmap(
    **plot_cache("quantms_clustered_heatmap", heatmap_z),
    id_col=id_col,
    data=heatmap_lazy,
    metadata=metadata_pl,
    row_cluster=True,
    col_cluster=True,
    title="Protein Abundance Heatmap (Z-score, clustered)",
    x_label="Samples",
    y_label="Proteins",
    colorscale=[[0, "#6699E0"], [0.5, "#FFFFFF"], [1, "#E06666"]],
    reversescale=False,
    intensity_label="Z-score",
    group_colors=group_colors,
)

state_manager = st.session_state.get("state")
# Scale height with the number of proteins so row labels stay readable -
# BaseComponent otherwise defaults to a flat 400px, too short for a
# dendrogram+heatmap composite with more than a handful of rows.
heatmap_height = max(600, min(1400, 300 + top_n * 20))
heatmap_component(state_manager=state_manager, height=heatmap_height)
st.caption(
    "Same Z-scored proteins as the heatmap, with rows and columns reordered so similar ones sit together "
    "(dendrograms). If the sample dendrogram splits along the group color bar, the groups differ consistently."
)

show_comparison_settings()

next_steps(
    ('enrichment', 'Find the biology behind groups of changed proteins.'),
    ('pca', 'Sample similarity at a glance.'),
    ('heatmap', 'The same proteins without clustering.'),
)
