import streamlit as st
from src.common.common import page_setup
from src.common.guide import page_intro
from src.WorkflowTest import WorkflowTest

params = page_setup()

wf = WorkflowTest()

st.title("Run the analysis")
page_intro(
    "Start the analysis and follow its progress. It runs on the server, so you can close "
    "this page and come back later with the same workspace. When it has finished, open "
    "**Run summary** under *Main results* in the sidebar.\n\n"
    "Each run is kept: to try other settings, change them and start again. The results "
    "pages show the newest run; switch runs with **Showing run** in the sidebar."
)
wf.show_execution_section()
