"""Run Overview Page: key figures and a per-file summary of the last run."""
from pathlib import Path

import streamlit as st

from src.common.common import page_setup
from src.common.insight_tables import show_insight_table
from src.common.overview import run_overview
from src.common.results_helpers import get_abundance_data, get_workflow_dir
from src.workflow.ParameterManager import ParameterManager

params = page_setup()
st.title("Run Overview")

st.markdown(
    """
Key figures of the last run and how each file did: MS2 spectra, identified PSMs
at the PSM FDR, peptides, precursor mass error, rescoring and search time.
"""
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

workflow_dir = get_workflow_dir(st.session_state["workspace"])
filter_dir = workflow_dir / "results" / "psm_filter"
cache_dir = workflow_dir / "results" / "insight_cache"
workflow_params = ParameterManager(workflow_dir, "TOPP Workflow").get_parameters_from_json()
mzml_files = [str(f) for f in workflow_params.get("mzML-files", [])]

if not filter_dir.is_dir() or not any(filter_dir.glob("*.idXML")):
    st.info("No results available yet. Please run the workflow first.")
    st.page_link("content/workflow_run.py", label="Go to Run", icon="🚀")
    st.stop()

log_file = workflow_dir / "logs" / "commands-and-run-times.log"
sources = sorted(filter_dir.glob("*.idXML")) + [log_file]


def _mtimes(paths) -> tuple:
    return tuple(Path(p).stat().st_mtime if Path(p).exists() else 0.0 for p in paths)


@st.cache_data(show_spinner="Summarizing the run...")
def _overview(workflow_dir: str, files: tuple, mtimes: tuple):
    """run_overview, cached on the inputs' and results' modification times."""
    return run_overview(workflow_dir, list(files))


per_file, figures = _overview(
    str(workflow_dir), tuple(mzml_files), _mtimes(mzml_files + [str(p) for p in sources])
)

abundance = get_abundance_data(st.session_state["workspace"])
proteins = samples = conditions = None
if abundance is not None:
    pivot_df, expr_df, group_map = abundance
    proteins = len(pivot_df)
    samples = len(expr_df.columns)
    conditions = len({str(g) for g in group_map.values() if str(g).strip()}) or None

fdr = float(workflow_params.get("psm-fdr-percent", 1.0))


def _show(value, unit: str = "") -> str:
    return "–" if value is None else f"{value:,}{unit}"


st.subheader("Key figures")
row1 = st.columns(4)
row1[0].metric("Quantified proteins", _show(proteins), help="Proteins in the quantification result.")
row1[1].metric("Peptides", _show(figures["Peptides"]), help=f"Distinct target peptide sequences at {fdr:g} % PSM FDR, over all files.")
row1[2].metric("Target PSMs", _show(figures["Target PSMs"]), help=f"Target PSMs at {fdr:g} % PSM FDR, over all files.")
row1[3].metric("PSMs per MS2 spectrum", _show(figures["PSMs per MS2 (%)"], " %"), help="Target PSMs divided by the MS2 spectra of all files.")
row2 = st.columns(4)
row2[0].metric(
    "Median precursor error",
    _show(figures["Median precursor error (ppm)"], " ppm"),
    help="Median m/z error of the target PSMs, after correcting picks of a 13C isotope peak.",
)
row2[1].metric("Files", _show(figures["Files"]))
row2[2].metric("Samples", _show(samples), help="Sample columns of the quantification result.")
row2[3].metric("Conditions", _show(conditions), help="Distinct groups set on the Configure page.")

st.subheader("Files")
st.caption(
    f"PSM counts are at {fdr:g} % PSM FDR. ID rate is target PSMs per MS2 spectrum. Mass error is "
    "the median precursor m/z error of the file's target PSMs. Percolator is whether it rescored "
    "the file. Search is CometAdapter's run time in seconds in the last run."
)
# One row per file; the fractional columns go as text, since the table stores
# floats as float32 and would show 36.8 as 36.79999923706055
per_file_view = per_file.copy()
for col in ("ID rate (%)", "Mass error (ppm)"):
    per_file_view[col] = per_file_view[col].map(lambda v: "" if v is None or v != v else f"{v:g}")
show_insight_table(
    per_file_view,
    name="overview_files",
    cache_dir=cache_dir,
    source_files=sources,
    title="run_overview",
    height=300,
    # Share the page width between the columns rather than folding the last ones away
    layout="fitColumns",
)

st.markdown("---")
st.markdown("**Next step:** Look at the identifications")
st.page_link("content/results_filtered.py", label="Go to Filtered PSMs", icon="🎯")
