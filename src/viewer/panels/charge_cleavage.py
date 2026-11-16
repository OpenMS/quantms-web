"""Charge states and missed cleavages as a Viewer component: a table of every
file and the share of each charge state and missed cleavage count per file."""

import hashlib
import shutil
from pathlib import Path

import pandas as pd
import polars as pl
import streamlit as st
from openms_insight import Table

from src.common.charge_cleavage import read_psms, stacked_shares, summary
from src.common.results_helpers import get_run_dir, get_run_params
from src.viewer.scope import needs_results, scoped, shared_state, show_fig


def _title(column: str) -> str:
    """A short column title: "2+" for charge 2, "≥1 MC" for one or more missed cleavages."""
    if column.startswith("Charge "):
        charge = column.split()[1]
        return f"≥{charge}" if charge.endswith("+") else f"{charge}+"
    if column.endswith(" missed (%)"):
        count = column.split()[0]
        return f"≥{count.rstrip('+')} MC" if count.endswith("+") else f"{count} MC"
    return {"Target PSMs": "PSMs", "Mean missed": "Mean MC"}.get(column, column)


def _column(column: str) -> dict:
    """Insight column definition: the short title, the full name as its tooltip, one decimal."""
    definition = {"field": column, "title": _title(column), "sorter": "number", "hozAlign": "right"}
    if column == "Target PSMs":
        return definition
    return {
        **definition,
        "formatter": "money",
        "formatterParams": {"precision": 2 if column == "Mean missed" else 1, "thousand": False},
        "headerTooltip": "Mean missed cleavages per peptide" if column == "Mean missed" else column,
    }


def render(ctx) -> None:
    run_dir = get_run_dir(ctx.workspace)
    filter_dir = run_dir / "results" / "psm_filter"
    cache_dir = run_dir / "results" / "insight_cache"
    idxmls = sorted(filter_dir.glob("*.idXML")) if filter_dir.is_dir() else []
    if not idxmls:
        return needs_results()
    mzml_names = tuple(Path(f).name for f in get_run_params(ctx.workspace).get("mzML-files", []) or [])
    mtimes = tuple(p.stat().st_mtime for p in idxmls)

    @st.cache_data(show_spinner="Counting charge states and missed cleavages...")
    def _summary(paths: tuple, names: tuple, mtimes: tuple):
        """summary() of the filtered PSMs, cached on the idXMLs' modification times."""
        read = [read_psms(Path(p), list(names)) for p in paths]
        psms = pd.concat([r[0] for r in read], ignore_index=True)
        enzymes = sorted({r[1] for r in read})
        return (*summary(psms), enzymes) if len(psms) else None

    result = _summary(tuple(str(p) for p in idxmls), mzml_names, mtimes)
    if result is None:
        return needs_results("No target PSMs in the filtered results.")
    per_file, charges, missed, enzymes = result
    st.caption(
        "Target PSMs at the PSM FDR. Most precursors of a tryptic digest carry two or three "
        "charges; many singly charged ones point to contamination. Missed cleavages are counted "
        f"per peptide for {', '.join(enzymes)}: a large share with one or more means an "
        "incomplete digestion. Shares in %; MC: missed cleavages."
    )

    key = hashlib.sha1("|".join(f"{p}:{m}" for p, m in zip(idxmls, mtimes)).encode()).hexdigest()[:10]
    cache_id = f"charge_cleavage_files_{key}"
    if not (cache_dir / cache_id).is_dir():
        for stale in cache_dir.glob("charge_cleavage_*"):
            shutil.rmtree(stale, ignore_errors=True)
        cache_dir.mkdir(parents=True, exist_ok=True)
        Table(
            cache_id=cache_id,
            data=pl.from_pandas(per_file).lazy(),
            cache_path=str(cache_dir),
            index_field="file_id",
            title="charge_and_missed_cleavages",
            layout="fitColumns",
            pagination=False,
            column_definitions=[{"field": "File", "title": "File", "widthGrow": 4}]
            + [_column(c) for c in per_file.columns if c not in ("file_id", "File")],
        )
    Table(cache_id=cache_id, cache_path=str(cache_dir))(
        key=scoped("charge_cleavage_files"), state_manager=shared_state(), height=90 + 34 * min(len(per_file), 12)
    )
    show_fig(
        stacked_shares(charges, per_file, "Precursor charge states of target PSMs", "Charge", "psms"),
        "qc-charge-states",
    )
    if missed.empty:
        st.caption("The search used no enzyme, so there are no missed cleavages to count.")
    else:
        show_fig(
            stacked_shares(missed, per_file, "Missed cleavages of identified peptides", "Missed cleavages", "peptides"),
            "qc-missed-cleavages",
        )
