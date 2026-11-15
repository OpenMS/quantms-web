"""GO enrichment of the proteins picked in the volcano plot or the group comparison table.

While nothing is picked (or with **All proteins**) it shows the GO annotations of
every tested protein: how many carry each term. Drawn by the template's GO panel
(``src/viz/go_panel.py``); this module is its source.
"""
import hashlib
from pathlib import Path

import streamlit as st

from src.common.differential import lead_accession
from src.common.guide import ensure_group_comparison
from src.common.results_helpers import get_abundance_data, get_id_column
from src.viewer import go_terms
from src.viewer.scope import GO_SET, heading, needs_results, scoped, shared_state
from src.viz.go_panel import show_go_panel

# Results kept per session (one per selection)
_RESULTS_KEPT = 8

HINT = (
    "Select proteins to test which GO terms are over-represented in them: **Select** in the toolbar of the "
    "volcano plot (All changed, Higher, Lower), box or lasso select on it, or tick rows of the group comparison table."
)


def _digest(*parts) -> str:
    return hashlib.sha1(repr(parts).encode()).hexdigest()[:12]


def _cached(digest: str, compute):
    store = st.session_state.setdefault("go_results", {})
    if digest not in store:
        store[digest] = compute()
        while len(store) > _RESULTS_KEPT:
            store.pop(next(iter(store)))
    return store[digest]


def render(ctx) -> None:
    abundance = get_abundance_data(ctx.workspace)
    if abundance is None:
        return needs_results()
    statistics_df = ensure_group_comparison(ctx.params)
    if statistics_df is None:
        return
    id_col = get_id_column(ctx.workspace, abundance[0])
    background = sorted({lead_accession(p) for p in statistics_df[id_col].dropna()})
    if not background:
        st.info("No protein was tested, so there are no GO terms to test.")
        return
    labels = {go_terms.clean(a): go_terms.label(a) for a in background}
    state = shared_state()
    selected = sorted({a for a in (state.get_selection(GO_SET) or []) if a})
    folder = Path(ctx.workspace, "insight_tables", "go")
    help_text = heading().get("help")

    def index():
        # The workflow looks the GO terms up after its run; only results of an
        # earlier build look them up here
        if go_terms.index_ready(folder, background):
            return go_terms.term_index(folder, background)
        with st.spinner(f"Looking up the GO terms of {len(background)} proteins at MyGene.info…"):
            return go_terms.term_index(folder, background)

    def enrichment():
        terms, error = index()
        if terms is None:
            return {}, error
        result = _cached(_digest(selected, background), lambda: go_terms.enrichment(terms, selected, background))
        return (result["categories"] if result["fg_count"] else None), error

    def counts():
        terms, error = index()
        if terms is None:
            return 0, {}, error
        annotated, tables = _cached(_digest("all", background), lambda: go_terms.term_counts(terms, background))
        return annotated, tables, error

    show_go_panel(
        scoped("go"), Path(ctx.workspace, "insight_tables"), state, GO_SET, bool(selected),
        f"**{len(selected)} of {len(background)} tested proteins selected.**", enrichment, counts, labels,
        hint=HINT, help=help_text,
    )
