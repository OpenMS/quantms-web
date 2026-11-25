"""The History of a run: on Explore results (the **History** button of the
experiment row) and Export results (a section of each run's card).

It shows how the run's experiment came about (``history.timeline``), the tool
settings that differ from the defaults, the processing report, the provenance
and **Start again with these settings**, which opens Describe experiment with
the experiment's files and the run's settings.
"""
from pathlib import Path

import pandas as pd
import streamlit as st

from src.common.provenance_ui import has_provenance, provenance
from src.viewer import history
from src.workflow.report_ui import report_buttons
from src.workflow.RunManager import RunManager

CONFIGURE_PAGE = Path("content", "workflow_configure.py")


def start_again(runs: RunManager, run_id: str) -> None:
    """Open Describe experiment with the settings ``run_id``'s experiment started
    with and all of its files."""
    from src.WorkflowTest import WorkflowTest  # the page set up the workspace already
    root = runs.experiment_root(run_id) or run_id
    files = runs.experiment_files(root, history.FILES_KEY)
    WorkflowTest().follow_up(root, changes={history.FILES_KEY: files} if files else None,
                             page=CONFIGURE_PAGE, note="to start again")


def history_section(runs: RunManager, run_id: str, key: str) -> None:
    """How ``run_id``'s experiment came about, its settings, report and provenance."""
    st.markdown("**How this result came about**")
    events = history.timeline(runs, run_id)
    st.markdown("\n".join(f"{i}. {history.event_markdown(e)}" for i, e in enumerate(events, start=1)))

    st.markdown("**Settings that differ from the defaults**")
    rows = history.changed_settings(runs, run_id)
    if rows:
        st.dataframe(pd.DataFrame(rows, columns=["Tool", "Parameter", "Default", "This result"]),
                     hide_index=True, width="stretch")
    else:
        st.caption("Every tool ran with its default settings.")

    members = runs.experiment(run_id) or [run_id]
    for member in members:
        if len(members) > 1:
            st.caption((runs.get(member) or {}).get("name", member))
        report_buttons(runs.workflow_dir, member, key=f"{key}-report-{member}")
    with st.container(horizontal=True, vertical_alignment="center"):
        if st.button("Start again with these settings", icon=":material/replay:", type="primary",
                     key=f"{key}-again",
                     help="Open Describe experiment with these files and settings. "
                          "The new run is kept next to this one."):
            start_again(runs, run_id)
    if has_provenance(runs, run_id):
        # Drawn inside an expander on Export results, so no expander of its own
        st.markdown("**Provenance**", help="The tools, versions and files of this run, as written to its "
                                           "QPX export.")
        provenance(runs, run_id)


OPEN_KEY = "viewer_history_open"


def _closed() -> None:
    st.session_state[OPEN_KEY] = False


@st.dialog("History", width="large", on_dismiss=_closed)
def history_dialog(runs: RunManager, run_id: str) -> None:
    st.markdown(f"#### {(runs.get(run_id) or {}).get('name', run_id)}")
    history_section(runs, run_id, key=f"history-dialog-{run_id}")
    st.caption("The history of every run is also on the card of the run on **Export results**.")


def history_button(runs: RunManager, run_id: str) -> None:
    """**History** of the run shown, for the experiment row of Explore results.

    The dialog stays open over reruns it did not ask for (its provenance table
    reruns the page when it first draws), until it is closed.
    """
    if st.button("History", icon=":material/history:", key="viewer_history",
                 help="How this result came about, its settings and its processing report"):
        st.session_state[OPEN_KEY] = True
    if st.session_state.get(OPEN_KEY):
        history_dialog(runs, run_id)
