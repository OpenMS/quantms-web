"""results_mass_shifts as a Viewer component."""

import hashlib
import shutil
from pathlib import Path
import polars as pl
import streamlit as st
from openms_insight import LinePlot, Table
from src.common.mass_shifts import MIN_PSMS, mass_shift_summary, psm_mass_shifts
from src.common.results_helpers import get_run_dir, get_run_params, search_tolerances
from src.viewer.scope import heading, needs_results, retitle, scoped, shared_state, zoom


def render(ctx) -> None:
    run_dir = get_run_dir(ctx.workspace)
    filter_dir = run_dir / "results" / "psm_filter"
    cache_dir = run_dir / "results" / "insight_cache"
    workflow_params = get_run_params(ctx.workspace)
    idxmls = sorted(filter_dir.glob("*.idXML")) if filter_dir.is_dir() else []
    if not idxmls:
        return needs_results()
    mzml_names = [Path(f).name for f in workflow_params.get("mzML-files", [])]
    mtimes = tuple(p.stat().st_mtime for p in idxmls)

    @st.cache_data(show_spinner="Computing mass shifts...")
    def _summary(paths: tuple, names: tuple, mtimes: tuple):
        """mass_shift_summary of the filtered PSMs, cached on the idXMLs' modification times."""
        return mass_shift_summary(psm_mass_shifts(list(paths), list(names)))

    shifts, histogram, psms = _summary(tuple(str(p) for p in idxmls), tuple(mzml_names), mtimes)
    tolerances = search_tolerances(workflow_params)
    tolerance, units = tolerances["precursor"], tolerances["precursor_unit"]
    if units == "Da" and tolerance and float(tolerance) >= 100:
        search_note = f"Open search: precursor tolerance ±{float(tolerance):g} Da."
    else:
        search_note = (
            "Closed search: only shifts within the precursor tolerance. "
            "A wider precursor tolerance in Configure finds modifications outside it."
        )
    if psms.empty:
        st.info("No target PSMs in the filtered results.")
        return
    modified = shifts[
        ~shifts["Interpretation"].isin(["Unmodified"]) & ~shifts["Interpretation"].str.contains("isotope")
    ]
    def _fixed(precision: int) -> dict:
        """Column options showing a number with ``precision`` decimals (the table stores float32)."""
        return {
            "sorter": "number",
            "hozAlign": "right",
            "formatter": "money",
            "formatterParams": {"precision": precision, "thousand": False},
        }

    # One set of caches per result: keyed on the idXMLs, their modification times
    # and the detected shifts (so a change in how shifts are found rebuilds them)
    key = hashlib.sha1(("|".join(f"{p}:{m}" for p, m in zip(idxmls, mtimes)) + shifts.to_csv()).encode()).hexdigest()[
        :10
    ]
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
            title="Shifts",
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
            title="PSMs of the selected shift",
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
    state_manager = shared_state()
    head = heading("Precursor mass shifts", note=search_note)
    retitle(LinePlot(cache_id=ids["histogram"], cache_path=str(cache_dir), **zoom(LinePlot, "mass_shifts")), head.get("title"))(
        key=scoped("mass_shift_histogram"), state_manager=state_manager, height=400, help=head.get("help")
    )
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
    retitle(Table(cache_id=ids["shifts"], cache_path=str(cache_dir), collapse_columns=False), "Shifts")(
        key=scoped("mass_shift_table"), state_manager=state_manager, height=140 + 34 * min(len(shifts), 10)
    )
    retitle(Table(cache_id=ids["psms"], cache_path=str(cache_dir), collapse_columns=False), "PSMs of the selected shift")(
        key=scoped("mass_shift_psm_table"), state_manager=state_manager, height=420
    )
