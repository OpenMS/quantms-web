"""Volcano Plot Results Page."""
import streamlit as st
import polars as pl
from src.common.common import page_setup
from src.common.differential import comparison_groups
from src.common.guide import ensure_group_comparison, needs_results, next_steps, page_intro, show_comparison_settings
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map, run_picker
from src.common.insight_tables import plot_cache
from openms_insight import VolcanoPlot

params = page_setup()
run_picker()
st.title("Volcano plot")

page_intro(
    "Every dot is a protein. The further **right or left**, the larger the change between your "
    "groups; the **higher**, the more confident the test is that the change is real. Proteins in "
    "the upper corners change strongly *and* reliably: these are your candidates.",
    background="""
- **x-axis, log2 fold change**: +1 means twice as abundant in the compared group, +2 four times,
  −1 half as abundant.
- **y-axis, −log10 adjusted p-value**: 1.3 corresponds to an adjusted p-value of 0.05, 2 to 0.01.
  The adjustment (Benjamini-Hochberg) accounts for testing thousands of proteins at once.
- **Cut-offs**: a protein is called up (red) or down (blue) when it passes both sliders.
  |log2 fold change| ≥ 1 and adjusted p ≤ 0.05 are common choices; loosen them to explore,
  tighten them for a short, high-confidence list.
- Hover over a dot to see the protein.
""",
    background_title="How to read a volcano plot",
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

result = get_abundance_data(st.session_state["workspace"])
if result is None:
    needs_results()

# Computed with the recommended settings unless the Fine-tune pages already did
statistics_df = ensure_group_comparison(params)
if statistics_df is None:
    st.stop()

if statistics_df.empty:
    st.info("No data available for volcano plot.")
    st.stop()

pivot_df, expr_df, group_map = result
id_col = get_id_column(st.session_state["workspace"], pivot_df)
sample_group_map = get_sample_group_map(st.session_state["workspace"], pivot_df, group_map)
direction = comparison_groups(sample_group_map, [c for c in pivot_df.columns if c in sample_group_map])
if direction:
    st.markdown(f"Right: higher in **{direction[1]}**. Left: higher in **{direction[0]}**.")

# 2. Clean data and convert to Polars for component input
volcano_df = statistics_df.dropna(subset=["log2FC", "p-adj"]).copy()
volcano_pl_lazy = pl.from_pandas(volcano_df).lazy()

# 3. Configure UI sliders (changing thresholds does not invalidate cache)
fc_thresh = st.slider(
    "Minimum change (|log2 fold change|)",
    min_value=0.5,
    max_value=3.0,
    value=1.0,
    step=0.1,
    help="1 = at least twice (or half) as abundant.",
)

p_thresh = st.slider(
    "Maximum adjusted p-value",
    min_value=0.001,
    max_value=0.1,
    value=0.05,
    step=0.001,
    help="0.05 means about 5 % of the proteins called significant are expected to be false positives.",
)

# 4. Initialize the OpenMS-Insight VolcanoPlot component
volcano_plot_component = VolcanoPlot(
    **plot_cache("quantms_volcano_plot", volcano_df),
    data=volcano_pl_lazy,
    log2fc_column="log2FC",
    pvalue_column="p-adj",
    label_column=id_col,
    up_color="#E74C3C",
    down_color="#3498DB",
    ns_color="#95A5A6",
    show_threshold_lines=True,
    threshold_line_style="dash",
)

# 5. Render the component
state_manager = st.session_state.get("state")  # Inject the project state management object

volcano_plot_component(
    state_manager=state_manager,
    fc_threshold=fc_thresh,
    p_threshold=p_thresh,
    max_labels=10,  # Display labels for the top N significant proteins
    height=600,
)

up_count = ((volcano_df["p-adj"] <= p_thresh) & (volcano_df["log2FC"] >= fc_thresh)).sum()
down_count = ((volcano_df["p-adj"] <= p_thresh) & (volcano_df["log2FC"] <= -fc_thresh)).sum()
if direction:
    st.markdown(f"**Higher in {direction[1]}:** {up_count} | **Lower in {direction[1]}:** {down_count}")
else:
    st.markdown(f"**Up-regulated:** {up_count} | **Down-regulated:** {down_count}")
show_comparison_settings()

next_steps(
    ("enrichment", "Find biological processes over-represented among these proteins."),
    ("compare", "The full table of compared proteins, most significant first."),
    ("pca", "Check that samples of the same group look alike."),
)
