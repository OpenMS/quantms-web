import streamlit as st
from src.common.common import page_setup
from src.common.guide import latest_empty_run, next_steps, page_intro, show_empty_run
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

wf_params = wf.parameter_manager.get_parameters_from_json()
selected = wf_params.get("mzML-files") or []
selected = [selected] if isinstance(selected, str) else selected
mode = wf_params.get("analysis-mode", "LFQ")
what = "Label-free analysis" if mode != "TMT" else "Isobaric (TMT/iTRAQ) analysis"
st.info(
    f"**{what}** of {len(selected)} run{'s' if len(selected) != 1 else ''}. "
    "Change it on **Describe experiment**.",
    icon="🧭",
)

wf.show_execution_section()

next_steps(
    ("overview", "See whether the run went well: identifications, error rate and quality."),
    ("compare", "Find the proteins that differ between your groups."),
)

