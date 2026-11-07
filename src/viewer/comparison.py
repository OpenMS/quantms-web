"""Edit comparison: the four steps of a run's group comparison, in a dialog."""
import json

import streamlit as st

from src.common.common import save_params
from src.common.differential import STEPS
from src.viewer.catalog import COMPARISON_STEPS, ViewerContext, draw
from src.viewer.scope import experiment, scoped


def _start_over(ctx: ViewerContext) -> None:
    """Forget the steps set by hand and go back to the recommended settings."""
    with experiment(ctx.run_id, ctx.index):
        for key in [*STEPS, "normalized_transform", "comparison_settings"]:
            st.session_state.pop(key, None)
    with open("default-parameters.json", "r", encoding="utf-8") as f:
        defaults = json.load(f)
    for key, value in defaults.items():
        if key.startswith("postproc-"):
            ctx.params[key] = value
            st.session_state.pop(key, None)
    save_params(ctx.params)


def _open_key(index: int) -> str:
    return f"comparison_dialog_open_{index}"


def open_comparison(index: int = 0) -> None:
    """Open the Edit comparison dialog of experiment ``index`` (drawn by ``edit_comparison_button``)."""
    st.session_state[_open_key(index)] = True


def edit_comparison_button(ctx: ViewerContext, button: bool = True) -> None:
    """A button opening the comparison steps of the experiment's run.

    Without ``button`` only the dialog is drawn, opened by ``open_comparison``
    (the Viewer's view bar does that for the first run). The dialog stays open
    over reruns of the page it did not ask for (the Viewer's components can
    rerun it right after the click), until it is closed or **Done** is pressed.
    """
    open_key = _open_key(ctx.index)

    def _closed() -> None:
        st.session_state[open_key] = False

    @st.dialog("Edit comparison", width="large", on_dismiss=_closed)
    def _dialog() -> None:
        tabs = st.tabs([label for _, label in COMPARISON_STEPS])
        for tab, (module, _) in zip(tabs, COMPARISON_STEPS):
            with tab:
                draw(module, ctx)
        left, right = st.columns(2)
        if left.button(
            "Start over with the recommended settings", key=scoped(f"comparison_reset_{ctx.index}"),
            help="Discard the steps changed here and compare the groups again with the defaults.",
        ):
            _start_over(ctx)
            st.rerun()
        if right.button("Done", type="primary", key=f"comparison_done_{ctx.index}"):
            _closed()
            st.rerun()

    if button and st.button("Edit comparison", key=f"comparison_open_{ctx.index}", icon=":material/tune:"):
        st.session_state[open_key] = True
    if st.session_state.get(open_key):
        _dialog()
