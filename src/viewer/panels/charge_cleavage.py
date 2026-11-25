"""Charge states and missed cleavages as a Viewer component: the library's
``charge_cleavage`` (a table of every file and two bar plots) on the run's
filtered PSMs."""

from pathlib import Path

import pandas as pd
import streamlit as st
from openms_insight import SeriesPlot, Table

from src.common.results_helpers import get_run_dir, get_run_params
from src.viewer.charge_cleavage_source import read_psms
from src.viewer.scope import heading, needs_results, scoped, shared_state, zoom
from src.viz.identification import charge_cleavage
from src.viz.plumbing import cache_name, component


@st.cache_data(show_spinner="Counting charge states and missed cleavages...")
def _summary(paths: tuple, names: tuple, mtimes: tuple):
    """The library summary of the filtered PSMs and the enzymes searched with,
    cached on the idXMLs' modification times (``mtimes`` only key the cache)."""
    read = [read_psms(Path(p), list(names)) for p in paths]
    psms = pd.concat([r[0] for r in read], ignore_index=True)
    enzymes = sorted({r[1] for r in read})
    return (*charge_cleavage.summary(psms), enzymes) if len(psms) else None


def _bars(name: str, settings: dict, cache_path: Path, **title) -> None:
    """One of the component's bar plots; the first takes the Viewer's heading.
    Cached by what it draws, so a change of the library's settings rebuilds it."""
    frame = settings.pop("data")
    SeriesPlot(
        cache_id=cache_name(f"quantms_{name}", frame, sorted(settings.items(), key=repr)),
        cache_path=str(cache_path), **zoom(SeriesPlot, name),
        **(settings | heading(**title)), data=frame.lazy(),
    )(key=scoped(name), state_manager=shared_state(), height=340)


def render(ctx) -> None:
    run_dir = get_run_dir(ctx.workspace)
    filter_dir = run_dir / "results" / "psm_filter"
    idxmls = sorted(filter_dir.glob("*.idXML")) if filter_dir.is_dir() else []
    if not idxmls:
        return needs_results()
    mzml_names = tuple(Path(f).name for f in get_run_params(ctx.workspace).get("mzML-files", []) or [])
    result = _summary(tuple(str(p) for p in idxmls), mzml_names, tuple(p.stat().st_mtime for p in idxmls))
    if result is None:
        return needs_results("No target PSMs in the filtered results.")
    per_file, charges, missed, enzymes = result
    cache_dir = run_dir / "results" / "insight_cache"

    note = (
        f"Target PSMs at the PSM FDR; missed cleavages per peptide for {', '.join(enzymes)}. "
        "Shares in %, MC: missed cleavages."
    )
    _bars("charge_cleavage", charge_cleavage.charge_settings(charges), cache_dir,
          title="Precursor charge states of target PSMs", note=note)
    if missed.empty:
        st.caption("The search used no enzyme, so there are no missed cleavages to count.")
    else:
        _bars("charge_cleavage_missed", charge_cleavage.cleavage_settings(missed), cache_dir,
              title="Missed cleavages of identified peptides")
    settings = charge_cleavage.table_settings(per_file)
    data = settings.pop("data")
    table = component(
        Table, cache_dir, cache_name("charge_cleavage_files", per_file, repr(settings)),
        lambda: {"data": data, **settings},
    )
    table(key=scoped("charge_cleavage_files"), state_manager=shared_state(),
          height=charge_cleavage.table_height(per_file))
