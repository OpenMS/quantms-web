"""Edit comparison: the four steps of a run's group comparison, in a dialog."""
import json

import streamlit as st

from src.common.common import save_params
from src.layout_manager import comparison as layout_comparison
from src.common.differential import STEPS
from src.viewer.catalog import COMPARISON_STEPS, ViewerContext, draw
from src.viewer.scope import experiment


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


def _key(index: int) -> str:
    return f"viewer_{index}"


def open_comparison(index: int = 0) -> None:
    """Open the Edit comparison dialog of experiment ``index`` (drawn by ``edit_comparison_button``)."""
    layout_comparison.open_comparison(_key(index))


def edit_comparison_button(ctx: ViewerContext, button: bool = True) -> None:
    """A button opening the comparison steps of the experiment's run, in the
    template's Edit comparison dialog (one tab per step).

    Without ``button`` only the dialog is drawn, opened by ``open_comparison``.
    """
    layout_comparison.edit_comparison(
        _key(ctx.index),
        [(label, lambda module=module: draw(module, ctx)) for module, label in COMPARISON_STEPS],
        on_reset=lambda: _start_over(ctx),
        reset_help="Discard the steps changed here and compare the groups again with the defaults.",
        button=button,
    )
