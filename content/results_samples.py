"""Samples Page: PSMs, peptides and proteins per sample, grouped by condition."""
from pathlib import Path

import plotly.express as px
import streamlit as st

from src.common.common import page_setup, show_fig
from src.common.insight_tables import show_insight_table
from src.common.results_helpers import get_abundance_data, get_sample_group_map, get_workflow_dir
from src.common.samples import COUNT_COLUMNS, channel_samples, condition_summary, lfq_samples
from src.workflow.ParameterManager import ParameterManager

params = page_setup()
st.title("Samples")

st.markdown(
    """
How many PSMs, peptides and proteins each sample contributed, grouped by the
conditions set on the Configure page. A sample far below the others in its
condition points to a failed run or a loading problem.
"""
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

workflow_dir = get_workflow_dir(st.session_state["workspace"])
results_dir = workflow_dir / "results"
cache_dir = results_dir / "insight_cache"
workflow_params = ParameterManager(workflow_dir, "TOPP Workflow").get_parameters_from_json()
analysis_mode = workflow_params.get("analysis-mode", "LFQ")
fdr = float(workflow_params.get("psm-fdr-percent", 1.0))


def _mtimes(paths) -> tuple:
    return tuple(Path(p).stat().st_mtime if Path(p).exists() else 0.0 for p in paths)


@st.cache_data(show_spinner="Counting per sample...")
def _lfq_samples(workflow_dir: str, files: tuple, groups: tuple, mtimes: tuple):
    """lfq_samples, cached on the group assignments and the results' modification times."""
    return lfq_samples(workflow_dir, list(files), dict(groups))


if analysis_mode == "LFQ":
    filter_dir = results_dir / "psm_filter"
    if not filter_dir.is_dir() or not any(filter_dir.glob("*.idXML")):
        st.info("No results available yet. Please run the workflow first.")
        st.page_link("content/workflow_run.py", label="Go to Run", icon="🚀")
        st.stop()
    mzml_files = [str(f) for f in workflow_params.get("mzML-files", [])]
    group_map = {
        key[len("mzML-group-"):]: value
        for key, value in workflow_params.items()
        if key.startswith("mzML-group-") and value
    }
    sources = sorted(filter_dir.glob("*.idXML")) + sorted((results_dir / "quant_results").glob("*.csv"))
    samples = _lfq_samples(
        str(workflow_dir), tuple(mzml_files), tuple(sorted(group_map.items())), _mtimes(sources)
    )
    st.caption(
        f"PSMs, peptides and proteins are the targets identified at {fdr:g} % PSM FDR, from the "
        "top hit of each spectrum. Quantified proteins and peptides have an intensity above 0 in "
        "the sample's ProteomicsLFQ result."
    )
else:
    abundance = get_abundance_data(st.session_state["workspace"])
    if abundance is None:
        st.info("No quantification result available yet. Please run the workflow first.")
        st.page_link("content/workflow_run.py", label="Go to Run", icon="🚀")
        st.stop()
    pivot_df, _, group_map = abundance
    sample_group_map = get_sample_group_map(st.session_state["workspace"], pivot_df, group_map)
    # The report's first columns are the protein and its description; the channels follow
    sample_columns = [c for c in pivot_df.columns[4:] if c != "PeptideSequence"]
    samples = channel_samples(pivot_df, sample_columns, sample_group_map)
    sources = sorted((results_dir / "quant_results").glob("*.csv"))
    st.caption(
        "Each sample is an isobaric channel. Quantified proteins have an intensity above 0 in the "
        "channel. Identifications are made per file, not per channel: see the Overview page."
    )

if samples.empty:
    st.info("No samples in the last run.")
    st.stop()

def _table_height(rows: int) -> int:
    """Room for every row up to a page of 12, plus the header and pager."""
    return 140 + 34 * min(rows, 12)


counts = [c for c in COUNT_COLUMNS if samples[c].notna().any()]
group_key = "|".join(f"{s}={c}" for s, c in zip(samples["Sample"], samples["Condition"]))

st.subheader("Per sample")
metric = st.selectbox("Count", counts, key="samples-metric")
fig = px.bar(
    samples,
    x="Sample",
    y=metric,
    color="Condition",
    category_orders={"Sample": samples["Sample"].tolist()},
)
# One full-width bar per sample: px gives each condition its own offset group,
# which would leave room beside every bar for the other conditions
fig.update_traces(offsetgroup=None, alignmentgroup=None)
fig.update_layout(xaxis_title=None, legend_title_text="Condition", bargap=0.2)
# Slant the sample names once they no longer fit side by side
label_chars = sum(len(s) for s in samples["Sample"])
fig.update_xaxes(tickangle=-35 if label_chars > 90 else 0, automargin=True)
show_fig(fig, f"samples_{metric.lower().replace(' ', '_')}")

show_insight_table(
    samples[["Sample", "Condition"] + counts],
    name="samples_per_sample",
    cache_dir=cache_dir,
    source_files=sources,
    extra_key=group_key,
    title="samples",
    height=_table_height(len(samples)),
    layout="fitColumns",
)

st.subheader("Per condition")
summary = condition_summary(samples)
# Means go as text, since the table stores floats as float32 and would show 36.8 as 36.79999923706055
for col in summary.columns[2:]:
    summary[col] = summary[col].map(lambda v: "" if v != v else f"{v:g}")
show_insight_table(
    summary,
    name="samples_per_condition",
    cache_dir=cache_dir,
    source_files=sources,
    extra_key=group_key,
    title="conditions",
    height=_table_height(len(summary)),
    layout="fitColumns",
)
if (samples["Condition"] == "Unassigned").any():
    st.caption("Samples without a condition are shown as Unassigned. Set conditions on the Configure page.")
    st.page_link("content/workflow_configure.py", label="Go to Configure", icon="⚙️")
