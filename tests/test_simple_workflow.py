import pytest
from streamlit.testing.v1 import AppTest
import json

"""
Tests for the Digest page functionality.

These tests verify:
- Page can be launched without errors
- Session state initialization works correctly
"""

@pytest.fixture
def launch():
    """Launch the Digest page for testing."""
    app = AppTest.from_file("content/digest.py")
    with open("settings.json", "r") as f:
        app.session_state.settings = json.load(f)
    app.session_state.settings["test"] = True
    app.secrets["workspace"] = "test"
    app.run(timeout=30)
    return app

def test_page_loads(launch):
    """Ensure the digest page loads without errors."""
    assert not launch.exception
