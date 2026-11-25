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
# Set while a run of an experiment of several runs is shown on its own (Single run)
SINGLE_KEY = "_viewer_single"


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

# The proteins picked for the GO enrichment: the volcano plot (box, lasso, its
# All changed / Higher / Lower buttons) and the check boxes of the group
# comparison table set them, the GO enrichment component tests them
GO_SET = "go_proteins"


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


def zoom(cls, name: str) -> dict:
    """``zoom_identifier`` for an Insight plot, so a view can keep its zoom.

    The plot then keeps its zoom in the experiment's state as ``zoom@<name>``
    (name: the catalog component), which 📌 stores with the view and ↺ applies.
    Empty with an openms-insight release whose plot has no such argument.
    """
    import inspect

    from src.layout_manager import zoom_identifier

    if "zoom_identifier" not in inspect.signature(cls.__init__).parameters:
        return {}
    return {"zoom_identifier": zoom_identifier(name)}


def change_cutoffs(params: dict) -> tuple[float, float]:
    """The cut-offs of a changed protein (Edit comparison, step 4), as now set.

    The values set in this session win over the saved ones, so the plots
    follow the dialog's inputs at once.
    """
    from src.common.differential import CHANGE_KEYS, change_cutoffs as saved

    return saved({**params, **{k: st.session_state[k] for k in CHANGE_KEYS if k in st.session_state}})


def shown_run() -> str | None:
    """The run of the experiment being drawn, None outside the Viewer."""
    return st.session_state.get(RUN_KEY)


def shown_single() -> bool:
    """Whether the run being drawn is shown on its own (Single run), not with its experiment."""
    return bool(st.session_state.get(SINGLE_KEY))


@contextmanager
def experiment(run_id: str, index: int, single: bool = False):
    """Draw components for ``run_id`` as experiment ``index``.

    ``single``: the run is shown on its own (Single run), so what compares
    samples covers its files only, not every run of its experiment.

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
    # A run on its own and its whole experiment have comparisons of their own
    chain = f"{run_id}#single" if single else run_id
    for key, value in chains.get(chain, {}).items():
        st.session_state[key] = value
    st.session_state[RUN_KEY] = run_id
    st.session_state[SCOPE_KEY] = index
    st.session_state[SINGLE_KEY] = single
    try:
        yield
    finally:
        chains[chain] = {k: st.session_state.pop(k) for k in RUN_DERIVED_KEYS if k in st.session_state}
        for key, value in (chains.get(outer_run, {}) if outer_run is not None else outer).items():
            st.session_state[key] = value
        st.session_state.pop(RUN_KEY, None)
        st.session_state.pop(SCOPE_KEY, None)
        st.session_state.pop(SINGLE_KEY, None)


# The catalog component being drawn, until its main Insight view claims its
# title and help (heading())
COMPONENT_KEY = "_viewer_component"


def heading(title: str | None = None, detail: str | None = None, note: str | None = None) -> dict:
    """``title=`` and ``help=`` for the main Insight view of the component being drawn.

    The first view of a component claims the component's label as its title,
    so the title sits inside the view, and its explanation as the view's
    info icon. Later views of the component (and views outside
    the Viewer) get ``title`` only. ``detail`` follows the label (``Volcano
    plot: B vs. A``); ``note`` (figures of the shown run) follows the
    explanation in the info icon.
    """
    name = st.session_state.pop(COMPONENT_KEY, None)
    if not name:
        return {"title": title} if title else {}
    from src.viewer.catalog import LABELS
    from src.viewer.teaching import help_text

    out = {"title": LABELS.get(name, title)}
    if detail and out["title"]:
        out["title"] = f"{out['title']}: {detail}"
    text = help_text(name)
    if text and note:
        text = f"{text}\n\n{note}"
    if text:
        out["help"] = text
    return out


def retitle(component, title: str | None):
    """Show a component rebuilt from its cache under ``title``."""
    if title and hasattr(component, "_title"):
        component._title = title
    return component


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
