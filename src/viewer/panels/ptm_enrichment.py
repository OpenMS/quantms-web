"""PTM enrichment: which protein modifications are enriched among the proteins
that change between the groups (src/psea).

The run computes it on the group comparison with the recommended settings when
PTM enrichment is on in Your analysis; changed in Edit comparison, or for a run
without it, the comparison is ranked again here and computed on request.
"""
from pathlib import Path

import pandas as pd
import streamlit as st

from src import psea
from src.common.guide import ensure_group_comparison
from src.common.insight_tables import show_insight_table
from src.common.results_helpers import get_abundance_data, get_id_column, get_results_dir
from src.viewer.scope import heading, needs_results, scoped, shown_run


def _column(table: pd.DataFrame, *names: str) -> str | None:
    """The first column whose name (case-insensitive) is one of ``names``, then one containing it."""
    lower = {c.lower(): c for c in table.columns}
    for name in names:
        if name.lower() in lower:
            return lower[name.lower()]
    for name in names:
        hit = next((c for c in table.columns if name.lower() in c.lower()), None)
        if hit:
            return hit
    return None


def _chart(table: pd.DataFrame, key: str) -> None:
    """The most enriched modifications as bars, coloured by their adjusted p-value."""
    import plotly.express as px

    from src.viewer.scope import show_fig

    term = _column(table, "PTM", "Term", "pathway")
    score = _column(table, "NES", "ES")
    padj = _column(table, "pvaladj", "padj", "p.adjust", "qvalue", "adj")
    if not term or not score or not pd.api.types.is_numeric_dtype(table[score]):
        return
    top = table.reindex(table[score].abs().sort_values(ascending=False).index).head(20)
    fig = px.bar(top.iloc[::-1], x=score, y=term, orientation="h",
                 color=padj if padj and pd.api.types.is_numeric_dtype(table[padj]) else None,
                 color_continuous_scale="Blues_r", height=max(260, 24 * len(top) + 80))
    fig.update_layout(yaxis_title=None, margin=dict(l=10, r=10, t=10, b=10))
    show_fig(fig, f"ptm_enrichment_{key}")


def _run_result(ctx, table: pd.DataFrame, organism: str, permutations: int) -> Path | None:
    """The run's own result, when it was computed on this ranking."""
    run = shown_run() or ctx.run_id
    if not run:
        return None
    folder = psea.run_results(get_results_dir(ctx.workspace, run))
    if folder is None:
        return None
    try:
        used = pd.read_csv(folder / psea.SCORES_FILE)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return None
    return folder if psea.digest(used, organism, permutations) == psea.digest(table, organism, permutations) else None


def render(ctx) -> None:
    abundance = get_abundance_data(ctx.workspace)
    if abundance is None:
        return needs_results()
    statistics_df = ensure_group_comparison(ctx.params)
    if statistics_df is None:
        return
    table = psea.scores(statistics_df, get_id_column(ctx.workspace, abundance[0]))
    if table.empty:
        st.info("No compared protein has a fold change, so there is nothing to rank.")
        return
    head = heading()
    st.markdown(f"**{head.get('title') or 'PTM enrichment'}**", help=head.get("help"))
    on, organism, permutations = psea.settings(ctx.params)
    if not ctx.params.get("ptm-organism"):
        organism = psea.guess_organism(statistics_df.iloc[:, 0].astype(str)) or organism

    c1, c2, c3 = st.columns([3, 1, 1], vertical_alignment="bottom")
    organism = c1.selectbox("Organism", psea.ORGANISMS, index=psea.ORGANISMS.index(organism)
                            if organism in psea.ORGANISMS else 0, key=scoped("ptm_organism"),
                            help=head.get("help"))
    permutations = int(c2.number_input("Permutations", 10, 100000, permutations, 100,
                                       key=scoped("ptm_permutations")))
    folder = _run_result(ctx, table, organism, permutations)
    cache = Path(ctx.workspace, "insight_tables", "ptm_enrichment")
    cached = cache / psea.digest(table, organism, permutations)
    if folder is None and (cached / psea.RESULTS_FILE).is_file():
        folder = cached
    if folder is None and c3.button("Compute", key=scoped("ptm_compute"), type="primary",
                                    use_container_width=True):
        with st.spinner(f"Testing the modifications of {len(table)} ranked proteins…"):
            folder, error = psea.compute(table, organism, permutations, cache)
        if error:
            st.session_state[scoped("ptm_error")] = error
    results = psea.read_results(folder)
    if results is None:
        error = st.session_state.get(scoped("ptm_error")) or psea.missing()
        if error:
            st.warning(f"PTM enrichment cannot be computed: {error}", icon="⚠️")
        st.caption(
            f"{len(table)} compared proteins, ranked by their log2 fold change. "
            + ("This run computed no PTM enrichment for this comparison: **Compute** tests it here. "
               if not on else "The comparison changed since the run: **Compute** tests it again. ")
            + "Switch **PTM enrichment** on in **Your analysis** on Describe experiment to have "
              "every run compute it."
        )
        return
    st.caption(f"{len(table)} compared proteins ranked by their log2 fold change, {organism}, "
               f"{permutations} permutations.")
    _chart(results, scoped("ptm"))
    show_insight_table(results, name=scoped("ptm_enrichment_table"), title="PTM enrichment",
                       height=360, extra_key=f"{folder}")
