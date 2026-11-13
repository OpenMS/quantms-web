import streamlit as st
from src.common.common import page_setup
from src.common.guide import PAGES, latest_empty_run, next_steps, page_intro, show_empty_run
from src.WorkflowTest import WorkflowTest
from src.workflow.RunManager import RUN_NAME_KEY, RUNNING

params = page_setup()

wf = WorkflowTest()

st.title("Run the analysis")
page_intro(
    "Start the analysis and follow its progress. It runs on the server, so you can close "
    "this page and come back later with the same workspace. When it has finished, press "
    "**Explore results**.\n\n"
    "Each run is kept: to try other settings, change them and start again. **Explore results** "
    "shows the newest run; pick another one in its run selector."
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
# The experiment's name: given beside Run analysis, or the run going now
running = wf.runs.list_runs(RUNNING)
name = wf_params.get(RUN_NAME_KEY) or (running[0].get("name", "") if running else "")
st.info(
    (f"**{name}**: " if name else "")
    + f"**{what}** of {len(selected)} run{'s' if len(selected) != 1 else ''}. "
    "Change it on **Describe experiment**.",
    icon="🧭",
)

wf.show_execution_section(results_page=PAGES["viewer"][0], results_label="Explore results")

next_steps(
    ("viewer", "See whether the run went well and which proteins differ between your groups."),
    ("export", "Keep the result files, the settings and the processing report of each run."),
)

