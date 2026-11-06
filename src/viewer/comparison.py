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


def edit_comparison_button(ctx: ViewerContext) -> None:
    """A button opening the comparison steps of the experiment's run."""

    @st.dialog("Edit comparison", width="large")
    def _dialog() -> None:
        st.caption(
            "The group comparison runs these four steps on the run's protein abundances, with "
            "recommended settings until you change one. Apply a step to update the steps after it, "
            "then press **Done** to see the result in the Viewer."
        )
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
            st.rerun()

    if st.button("Edit comparison", key=f"comparison_open_{ctx.index}", icon=":material/tune:"):
        _dialog()
