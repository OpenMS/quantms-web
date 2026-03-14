import pytest
from streamlit.testing.v1 import AppTest
import json

@pytest.fixture
def launch():
    """Launch the Workflow Configure Streamlit page for testing."""
    app = AppTest.from_file("content/workflow_configure.py")
    with open("settings.json", "r") as f:
        app.session_state.settings = json.load(f)
    app.session_state.settings["test"] = True
    app.secrets["workspace"] = "test"
    app.run(timeout=10)
    return app

def test_page_loads(launch):
    """Ensure the workflow configure page loads without errors."""
    assert not launch.exception