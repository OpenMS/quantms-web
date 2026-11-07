"""Export results (local app) or Download results (online): the files of every run."""
from pathlib import Path

import streamlit as st

from src.common import result_export as rx
from src.common.common import export_folder_input, page_setup
from src.common.guide import needs_results
from src.common.provenance_ui import has_provenance, provenance
from src.common.results_helpers import get_runs
from src.workflow.report_ui import report_buttons
from src.workflow.RunManager import FINISHED, MAX_NAME_LENGTH, RUNNING, format_time
from src.workflow.run_widgets import run_import

page_setup()
local = rx.is_local()
st.title("Export results" if local else "Download results")
st.markdown(
    "One folder per run: its results, settings and logs." if local
    else "One ZIP per run: its results, settings and logs."
)

runs = get_runs(Path(st.session_state["workspace"]))
all_runs = runs.list_runs()
finished = [r["id"] for r in all_runs if r.get("status") == FINISHED]
if not all_runs:
    run_import(runs)
    needs_results()

export_dir = export_folder_input("results-export") if local else None
if finished:
    rx.all_runs_button(runs, finished, export_dir)

for i, run in enumerate(all_runs):
    run_id = run["id"]
    name = run.get("name", run_id)
    status = run.get("status")
    files = rx.run_files(runs, run_id)
    editing = f"_export-rename-{run_id}"
    with st.container(border=True):
        name_col, button_col, delete_col = st.columns([3, 1.2, 0.5], vertical_alignment="center")
        with name_col:
            if st.session_state.get(editing):
                with st.container(horizontal=True, vertical_alignment="bottom"):
                    st.text_input("New name", value=name, max_chars=MAX_NAME_LENGTH,
                                  key=f"export-name-{run_id}")
                    if st.button("Save", key=f"export-save-{run_id}", icon=":material/save:"):
                        if runs.rename(run_id, st.session_state.get(f"export-name-{run_id}", "")):
                            st.session_state[editing] = False
                            st.rerun()
                        st.error("The name can not be empty.")
            else:
                with st.container(horizontal=True, vertical_alignment="center"):
                    st.markdown(f"**{name}**")
                    if st.button("", key=f"export-edit-{run_id}", icon=":material/edit:",
                                 help="Rename", type="tertiary"):
                        st.session_state[editing] = True
                        st.rerun()
            root = runs.experiment_root(run_id)
            adds = (f" · adds files to {(runs.get(root) or {}).get('name', root)}"
                    if root and root != run_id else "")
            state = "" if status == FINISHED else f" · {status}"
            st.caption(f"{format_time(run.get('created'))} · {len(files)} files, "
                       f"{rx.format_size(rx.total_size(files))}{adds}{state}")
        message_area = st.container()
        with button_col:
            if status == FINISHED:
                rx.run_button(runs, run_id, str(i), export_dir,
                              label=None if local else "Download ZIP", message_area=message_area)
        with delete_col:
            with st.popover("", icon=":material/delete:", help=f"Delete {name} and its files",
                            width="stretch", disabled=status == RUNNING):
                st.markdown(f"Delete **{name}** and all its files?")
                if st.button("Delete", key=f"export-delete-{run_id}", type="primary"):
                    runs.delete(run_id)
                    st.rerun()
        with st.expander("Files"):
            rx.file_table(files)
        if status == FINISHED:
            with st.expander("Processing report"):
                report_buttons(runs.workflow_dir, run_id, key=f"export-report-{run_id}")
            if has_provenance(runs, run_id):
                with st.expander("Provenance"):
                    provenance(runs, run_id)

run_import(runs)
