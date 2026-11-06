"""What one experiment of the Viewer shows: its run, its widget keys, its comparison.

The Viewer draws several experiments on one page, each showing a run of its own.
The results code reads the shown run through ``results_helpers`` and keeps the
group comparison (filtered, imputed, normalized and tested tables) in session
state, so ``experiment()`` points both at the experiment's run while its
components are drawn, and ``scoped()`` gives each experiment its own widget keys.
"""
from contextlib import contextmanager

import streamlit as st

from src.common.guide import page_link

# Session-state keys: the run shown, the experiment drawn, and the comparison of every run
RUN_KEY = "_viewer_run"
SCOPE_KEY = "_viewer_scope"
CHAINS_KEY = "_viewer_chains"


def scoped(key: str) -> str:
    """``key`` for the experiment being drawn and its run.

    The run is part of the key so that picking another run draws new components:
    an Insight table kept under one key keeps the columns of the run it showed
    first, and a selection made in one run means nothing in another.
    """
    index = st.session_state.get(SCOPE_KEY)
    run = st.session_state.get(RUN_KEY)
    if run is not None:
        key = f"{key}_{run}"
    return f"{key}_e{index}" if index else key


# Identifiers of the PSM components that are row numbers of their own table,
# so they mean another PSM in each (database search, rescored, filtered)
PSM_PRIVATE = ("identification", "peak")


def shared_state(owner: str = "", private: tuple = ()):
    """The Insight state every component of the experiment being drawn shares.

    One per experiment and run, not per component or view: a protein, peptide
    or spectrum picked in one component stays picked in the others, also after
    switching views. ``private`` identifiers (a row number of the component's
    own table) are kept apart for ``owner``.
    """
    from src.layout_manager import experiment_state

    index = st.session_state.get(SCOPE_KEY) or 0
    return experiment_state("viewer", index, st.session_state.get(RUN_KEY), private=private, owner=owner)


def shown_run() -> str | None:
    """The run of the experiment being drawn, None outside the Viewer."""
    return st.session_state.get(RUN_KEY)


@contextmanager
def experiment(run_id: str, index: int):
    """Draw components for ``run_id`` as experiment ``index``.

    The group comparison of the run is put in place of whatever comparison the
    session holds, and saved back for the run afterwards, so two experiments
    showing two runs never mix their comparisons.
    """
    from src.common.results_helpers import RUN_DERIVED_KEYS

    chains = st.session_state.setdefault(CHAINS_KEY, {})
    outer_run = st.session_state.get("_results_run")
    outer = {k: st.session_state.pop(k) for k in RUN_DERIVED_KEYS if k in st.session_state}
    if outer_run is not None:
        chains[outer_run] = outer
    for key, value in chains.get(run_id, {}).items():
        st.session_state[key] = value
    st.session_state[RUN_KEY] = run_id
    st.session_state[SCOPE_KEY] = index
    try:
        yield
    finally:
        chains[run_id] = {k: st.session_state.pop(k) for k in RUN_DERIVED_KEYS if k in st.session_state}
        for key, value in (chains.get(outer_run, {}) if outer_run is not None else outer).items():
            st.session_state[key] = value
        st.session_state.pop(RUN_KEY, None)
        st.session_state.pop(SCOPE_KEY, None)


def needs_results(
    message: str = "No results of this kind in this run.",
) -> None:
    """Say that the run has nothing to show here (the component then returns)."""
    st.info(message)


def show_fig(fig, download_name: str, container_width: bool = True, selection_session_state_key: str = "") -> None:
    """``common.show_fig`` with a key per experiment, so one plot can be shown twice."""
    from src.common.common import follow_app_theme, show_fig as common_show_fig

    if selection_session_state_key or not st.session_state.get(SCOPE_KEY):
        common_show_fig(fig, download_name, container_width, scoped(selection_session_state_key) if selection_session_state_key else "")
        return
    follow_app_theme(fig)
    st.plotly_chart(
        fig,
        use_container_width=container_width,
        key=scoped(f"fig_{download_name}"),
        config={
            "displaylogo": False,
            "toImageButtonOptions": {"filename": download_name, "format": st.session_state.get("image-format", "svg")},
        },
    )


__all__ = ["experiment", "needs_results", "page_link", "scoped", "shown_run", "show_fig"]
