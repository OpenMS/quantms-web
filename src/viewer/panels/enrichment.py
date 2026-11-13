"""GO enrichment of the proteins picked in the volcano plot or the group comparison table.

While nothing is picked it shows the GO annotations of every tested protein:
how many carry each term.
"""
import hashlib
from pathlib import Path

import streamlit as st

from src.common.differential import lead_accession
from src.common.guide import ensure_group_comparison, show_comparison_settings
from src.common.insight_tables import show_insight_table, show_series
from src.common.results_helpers import get_abundance_data, get_id_column
from src.viewer import go_terms
from src.viewer.scope import GO_SET, needs_results, scoped, shared_state

TERM_COUNTS = [10, 15, 20, 30, 50]
# Results kept per session (one per selection)
_RESULTS_KEPT = 8

HINT = (
    "Pick proteins to test which GO terms are over-represented in them: use **All changed**, "
    "**Higher** or **Lower** above the volcano plot, box or lasso select in its mode bar (Shift "
    "adds), or the check boxes of the group comparison table (its check mark button turns them on)."
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


def _clear(state_manager) -> None:
    # Set to None, not removed: the components keep a selection Python no
    # longer sends and would send it back
    state_manager.set_selection(GO_SET, None)


def _pick_terms():
    pick, size = st.columns([4, 1], vertical_alignment="bottom")
    with pick:
        category = st.segmented_control(
            "GO category", list(go_terms.CATEGORIES), default="BP", key=scoped("go_category"),
            format_func=go_terms.CATEGORIES.get, label_visibility="collapsed",
        ) or "BP"
    count = size.selectbox(
        "Terms in the chart", TERM_COUNTS, index=1, key=scoped("go_term_count"),
        help="How many terms the bar chart shows. The table lists all of them.",
    )
    return category, count


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

    if selected:
        head, clear = st.columns([5, 1], vertical_alignment="center")
        head.markdown(f"**{len(selected)} of {len(background)} tested proteins picked.**")
        clear.button("Clear selection", key=scoped("go_clear"), on_click=_clear, args=(state,), use_container_width=True)
    else:
        st.info(HINT)

    message = (
        "Testing GO terms…" if go_terms.index_ready(folder, background)
        else f"Looking up the GO terms of {len(background)} proteins at MyGene.info…"
    )
    with st.spinner(message):
        index, error = go_terms.term_index(folder, background)
    if error:
        st.warning(error + " Try again when the network is back.")
    if index is None:
        return

    if not selected:
        annotated, tables = _cached(_digest("all", background), lambda: go_terms.term_counts(index, background))
        if not annotated:
            st.info("None of the tested proteins has a GO annotation at MyGene.info.")
            return
        category, count = _pick_terms()
        table = tables[category]
        if not table.height:
            st.info(f"No tested protein has {go_terms.CATEGORIES[category].lower()} terms.")
            return
        show_series(
            go_terms.count_bars(table, count), name=f"go-counts-{category}",
            x_column="term", y_column="proteins", mode="bars", orientation="h", show_legend=False,
            title=f"{go_terms.CATEGORIES[category]}: the most frequent terms", x_label="GO term", y_label="Proteins",
            height=max(300, 60 + 24 * min(count, table.height)),
        )
        show_insight_table(
            go_terms.count_table(table, labels, annotated), name=f"go-counts-table-{category}",
            title=f"GO terms of all tested proteins ({table.height})", height=360, layout="fitColumns",
            column_overrides={"% of annotated": {"formatter": "money", "formatterParams": {"thousand": False, "precision": 1}}},
        )
        st.caption(
            f"GO annotations of all {len(background)} tested proteins: {annotated} have GO terms "
            "(MyGene.info, by UniProt accession)."
        )
        show_comparison_settings()
        return

    result = _cached(_digest(selected, background), lambda: go_terms.enrichment(index, selected, background))
    if not result["fg_count"]:
        st.info("None of the picked proteins has a GO annotation at MyGene.info.")
        return
    category, count = _pick_terms()
    table = result["categories"].get(category)
    if table is None or not len(table):
        st.info(f"The picked proteins have no {go_terms.CATEGORIES[category].lower()} terms.")
    else:
        bars = go_terms.term_bars(table, count)
        show_series(
            bars, name=f"go-enriched-{category}",
            x_column="term", y_column="score", mode="bars", orientation="h", show_legend=False,
            hover_column="hover", title=f"{go_terms.CATEGORIES[category]}: the most significant terms",
            x_label="GO term", y_label="-log10(adjusted p-value)", height=max(300, 60 + 24 * len(bars)),
        )
        terms = go_terms.term_table(table, labels)
        show_insight_table(
            terms, name=f"go-terms-{category}", title=f"GO terms ({len(terms)})", height=360, layout="fitColumns",
            column_overrides={"Fold": {"formatter": "money", "formatterParams": {"thousand": False, "precision": 2}}},
        )
    st.caption(
        f"{result['fg_count']} of {result['fg_total']} picked and {result['bg_count']} of {result['bg_total']} "
        "tested proteins have GO terms (MyGene.info, by UniProt accession). One-sided Fisher's exact test, "
        "Benjamini-Hochberg adjusted per category."
    )
    show_comparison_settings()
