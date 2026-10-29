"""Mass Shifts Page: precursor mass shifts of the identified PSMs and their Unimod matches."""
import hashlib
import shutil
from pathlib import Path

import polars as pl
import streamlit as st
from openms_insight import LinePlot, StateManager, Table

from src.common.common import page_setup
from src.common.mass_shifts import MIN_PSMS, mass_shift_summary, psm_mass_shifts
from src.common.results_helpers import get_workflow_dir
from src.workflow.ParameterManager import ParameterManager

params = page_setup()
st.title("Mass Shifts")

st.markdown(
    """
The difference between each identified PSM's observed precursor mass and the
theoretical mass of its peptide. With the **Open Search** preset (±500 Da
precursor tolerance) a peptide can be identified with a modification the search
did not include, and recurring shifts reveal which ones. Each shift is matched
to Unimod. Click a peak or a row to list its PSMs.
"""
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

workflow_dir = get_workflow_dir(st.session_state["workspace"])
filter_dir = workflow_dir / "results" / "psm_filter"
cache_dir = workflow_dir / "results" / "insight_cache"
workflow_params = ParameterManager(workflow_dir, "TOPP Workflow").get_parameters_from_json()

idxmls = sorted(filter_dir.glob("*.idXML")) if filter_dir.is_dir() else []
if not idxmls:
    st.info("No results available yet. Please run the workflow first.")
    st.page_link("content/workflow_run.py", label="Go to Run", icon="🚀")
    st.stop()

mzml_names = [Path(f).name for f in workflow_params.get("mzML-files", [])]
mtimes = tuple(p.stat().st_mtime for p in idxmls)


@st.cache_data(show_spinner="Computing mass shifts...")
def _summary(paths: tuple, names: tuple, mtimes: tuple):
    """mass_shift_summary of the filtered PSMs, cached on the idXMLs' modification times."""
    return mass_shift_summary(psm_mass_shifts(list(paths), list(names)))


shifts, histogram, psms = _summary(tuple(str(p) for p in idxmls), tuple(mzml_names), mtimes)

comet = workflow_params.get("CometAdapter", {})
tolerance = comet.get("precursor_mass_tolerance")
units = comet.get("precursor_error_units")
if units == "Da" and tolerance and float(tolerance) >= 100:
    st.caption(f"Open search: precursor tolerance ±{float(tolerance):g} Da.")
else:
    st.caption(
        "This run used a closed search, so only shifts within the precursor tolerance appear "
        "(0 Da and isotope picks). Apply the **Open Search** preset on the Configure page to "
        "look for unexpected modifications."
    )
    st.page_link("content/workflow_configure.py", label="Go to Configure", icon="⚙️")

if psms.empty:
    st.info("No target PSMs in the filtered results.")
    st.stop()

modified = shifts[~shifts["Interpretation"].isin(["Unmodified"]) & ~shifts["Interpretation"].str.contains("isotope")]
cols = st.columns(4)
cols[0].metric("Target PSMs", f"{len(psms):,}")
cols[1].metric("Shifts", f"{len(shifts):,}", help=f"Recurring shifts with at least {MIN_PSMS} PSMs.")
cols[2].metric(
    "PSMs with a modification shift",
    f"{int(modified['PSMs'].sum()):,}",
    help="PSMs in shifts other than 0 Da and isotope picks.",
)
cols[3].metric(
    "Unexplained shifts",
    f"{int((shifts['Interpretation'] == '').sum()):,}",
    help="Shifts with no Unimod modification within 0.02 Da.",
)


def _fixed(precision: int) -> dict:
    """Column options showing a number with ``precision`` decimals (the table stores float32)."""
    return {"sorter": "number", "hozAlign": "right", "formatter": "money",
            "formatterParams": {"precision": precision, "thousand": False}}


# One set of caches per result: keyed on the idXMLs, their modification times
# and the detected shifts (so a change in how shifts are found rebuilds them)
key = hashlib.sha1(
    ("|".join(f"{p}:{m}" for p, m in zip(idxmls, mtimes)) + shifts.to_csv()).encode()
).hexdigest()[:10]
ids = {name: f"mass_shift_{name}_{key}" for name in ("histogram", "shifts", "psms")}
if not all((cache_dir / cache_id).is_dir() for cache_id in ids.values()):
    for stale in cache_dir.glob("mass_shift_*"):
        shutil.rmtree(stale, ignore_errors=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    LinePlot(
        cache_id=ids["histogram"],
        data=pl.from_pandas(histogram).lazy(),
        cache_path=str(cache_dir),
        x_column="delta",
        y_column="psms",
        title="Precursor mass shifts",
        x_label="Mass shift (Da)",
        y_label="PSMs",
        highlight_column="peak",
        annotation_column="label",
        interactivity={"mass_shift": "shift_id"},
    )
    Table(
        cache_id=ids["shifts"],
        data=pl.from_pandas(shifts).lazy(),
        cache_path=str(cache_dir),
        interactivity={"mass_shift": "shift_id"},
        index_field="shift_id",
        title="mass_shifts",
        layout="fitColumns",
        column_definitions=[
            {"field": "Observed Δ (Da)", "title": "Observed Δ (Da)", **_fixed(4)},
            {"field": "PSMs", "title": "PSMs", "sorter": "number", "hozAlign": "right"},
            {"field": "Interpretation", "title": "Interpretation"},
            {"field": "Unimod Δ (Da)", "title": "Unimod Δ (Da)", **_fixed(4)},
            {"field": "Deviation (mDa)", "title": "Deviation (mDa)", **_fixed(1)},
            {"field": "Residues", "title": "Residues"},
            {"field": "Other matches", "title": "Other matches"},
        ],
    )
    psm_data = pl.from_pandas(psms.reset_index(drop=True)).with_row_index("psm_id")
    Table(
        cache_id=ids["psms"],
        data=psm_data.lazy(),
        cache_path=str(cache_dir),
        filters={"mass_shift": "shift_id"},
        index_field="psm_id",
        title="mass_shift_psms",
        initial_sort=[{"column": "score", "dir": "asc"}],
        column_definitions=[
            {"field": "sequence", "title": "Sequence"},
            {"field": "delta", "title": "Δ (Da)", **_fixed(4)},
            {"field": "charge", "title": "Z", "sorter": "number"},
            {"field": "mz", "title": "m/z", **_fixed(4)},
            {"field": "rt", "title": "RT", **_fixed(1)},
            {"field": "score", "title": "Score", "sorter": "number"},
            {"field": "file", "title": "File"},
        ],
    )

state_manager = StateManager()
LinePlot(cache_id=ids["histogram"], cache_path=str(cache_dir))(
    key="mass_shift_histogram", state_manager=state_manager, height=400
)

st.subheader("Shifts")
Table(cache_id=ids["shifts"], cache_path=str(cache_dir))(
    key="mass_shift_table", state_manager=state_manager, height=140 + 34 * min(len(shifts), 10)
)

st.subheader("PSMs of the selected shift")
Table(cache_id=ids["psms"], cache_path=str(cache_dir))(
    key="mass_shift_psm_table", state_manager=state_manager, height=420
)
