"""OpenMS-Insight image downloads follow the sidebar's "image export format".

page_setup() hands the selectbox value to openms_insight.set_image_format(),
which keeps it in the session's state (under the key below) for every
Insight component the page draws.
"""

import json
import sys
import types

import pytest
from streamlit.testing.v1 import AppTest

pytest.importorskip("openms_insight")
from openms_insight.core import image_format  # noqa: E402

PAGE = "content/workflow_run.py"


def _app() -> AppTest:
    app = AppTest.from_file(PAGE)
    with open("settings.json", "r", encoding="utf-8") as f:
        app.session_state.settings = json.load(f)
    app.session_state.settings["test"] = True
    app.secrets["workspace"] = "test"
    return app


def test_sidebar_format_reaches_insight():
    app = _app()
    app.run(timeout=60)
    assert not app.exception
    assert app.session_state[image_format._SESSION_KEY] == "svg"

    app.selectbox(key="image-format").set_value("png").run(timeout=60)
    assert not app.exception
    assert app.session_state[image_format._SESSION_KEY] == "png"


@pytest.fixture
def app():
    return _app()


def test_options_follow_insight_list(app, monkeypatch):
    # The sidebar offers exactly the formats Insight lists, in its order
    fake = types.ModuleType("openms_insight")
    fake.IMAGE_FORMATS = ("png", "svg")
    fake.set_image_format = lambda fmt: None
    monkeypatch.setitem(sys.modules, "openms_insight", fake)
    app.run(timeout=60)
    assert not app.exception
    box = app.selectbox(key="image-format")
    assert box.options == ["png", "svg"]
    assert box.value == "svg"


def test_saved_format_insight_dropped(app, monkeypatch):
    # A saved format Insight no longer lists falls back to its first one
    fake = types.ModuleType("openms_insight")
    fake.IMAGE_FORMATS = ("png",)
    fake.set_image_format = lambda fmt: None
    monkeypatch.setitem(sys.modules, "openms_insight", fake)
    app.run(timeout=60)
    assert not app.exception
    assert app.selectbox(key="image-format").value == "png"
