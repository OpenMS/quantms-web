"""Run Overview Page: key figures and a per-file summary of the last run."""
from pathlib import Path

import streamlit as st

from src.common.common import page_setup
from src.common.guide import needs_results, next_steps, page_intro
from src.common.insight_tables import show_insight_table
from src.common.overview import run_overview
from src.common.results_helpers import get_abundance_data, run_picker, get_run_dir, get_run_params

params = page_setup()
run_picker()
st.title("Run summary")

page_intro(
    "Did the run work? The key figures below say how many proteins and peptides were found, "
    "and the file table shows whether any sample stands out. Start here after every run.",
    background="""
- **PSM** (peptide-spectrum match): one fragmentation spectrum explained by one peptide.
  Every protein is backed by peptides, and every peptide by one or more PSMs.
- **PSMs per MS2 spectrum** (identification rate): the share of fragmentation spectra that
  were explained by a peptide. 20 to 50 % is typical for a good tryptic digest on a modern
  instrument; much lower values point to the wrong database, wrong modifications or wrong
  mass tolerances.
- **Precursor mass error**: how far measured peptide masses lie from the expected ones, in
  parts per million (ppm). On an Orbitrap the median should be within a few ppm of 0;
  a large offset means the instrument needs calibration.
- A **file far below the others** in PSMs or ID rate usually had a problem during sample
  preparation or measurement, and may distort the group comparison.
""",
    background_title="How to read these numbers",
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

run_dir = get_run_dir(st.session_state["workspace"])
filter_dir = run_dir / "results" / "psm_filter"
cache_dir = run_dir / "results" / "insight_cache"
workflow_params = get_run_params(st.session_state["workspace"])
mzml_files = [str(f) for f in workflow_params.get("mzML-files", [])]

if not filter_dir.is_dir() or not any(filter_dir.glob("*.idXML")):
    needs_results()

log_file = run_dir / "logs" / "commands-and-run-times.log"
sources = sorted(filter_dir.glob("*.idXML")) + [log_file]


def _mtimes(paths) -> tuple:
    return tuple(Path(p).stat().st_mtime if Path(p).exists() else 0.0 for p in paths)


@st.cache_data(show_spinner="Summarizing the run...")
def _overview(workflow_dir: str, files: tuple, mtimes: tuple):
    """run_overview, cached on the inputs' and results' modification times."""
    return run_overview(workflow_dir, list(files))


per_file, figures = _overview(
    str(run_dir), tuple(mzml_files), _mtimes(mzml_files + [str(p) for p in sources])
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
row2[3].metric("Conditions", _show(conditions), help="Distinct groups set on the Describe experiment page.")

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

next_steps(
    ("compare", "See which proteins differ between your groups."),
    ("proteins", "Look up a protein of interest and the peptides that support it."),
    ("download", "Get the protein and peptide tables for your own analysis."),
)
