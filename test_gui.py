import json

import pytest
from streamlit.testing.v1 import AppTest

# Pages that AppTest.from_file can load in isolation. Pages using st.page_link
# require streamlit's navigation context (only set up when app.py runs), so they
# are covered indirectly by test_app_loads below.
DIRECTLY_TESTABLE_PAGES = [
    "content/workflow_fileupload.py",
    "content/workflow_configure.py",
    "content/workflow_run.py",
    "content/results_library.py",
    "content/results_proteomicslfq.py",
]


def _init(apptest):
    with open("settings.json", "r") as f:
        apptest.session_state.settings = json.load(f)
    apptest.session_state.settings["test"] = True
    apptest.secrets["workspace"] = "test"
    return apptest


@pytest.fixture
def launch(request):
    return _init(AppTest.from_file(request.param))


@pytest.mark.parametrize("launch", DIRECTLY_TESTABLE_PAGES, indirect=True)
def test_page_loads(launch):
    launch.run(timeout=30)
    assert not launch.exception


def test_app_loads():
    app = _init(AppTest.from_file("app.py"))
    app.run(timeout=30)
    assert not app.exception
