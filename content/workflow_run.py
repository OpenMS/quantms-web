import streamlit as st
from src.common.common import page_setup
from src.common.guide import latest_empty_run, page_intro, show_empty_run
from src.WorkflowTest import WorkflowTest

params = page_setup()

wf = WorkflowTest()

st.title("Run the analysis")
page_intro(
    "Start the analysis and follow its progress. It runs on the server, so you can close "
    "this page and come back later with the same workspace. When it has finished, open "
    "**Run summary** under *Main results* in the sidebar.\n\n"
    "Each run is kept: to try other settings, change them and start again. The results "
    "pages show the newest run; switch runs with **Showing run** at the top of each results page."
)
# A run that kept no PSM or protein stops with an error; say why in plain words
if "workspace" in st.session_state:
    empty = latest_empty_run(st.session_state["workspace"])
    if empty is not None:
        _, run_name, diagnosis = empty
        show_empty_run(run_name, diagnosis, links=False)

wf.show_execution_section()

