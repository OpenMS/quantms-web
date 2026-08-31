"""Smoke tests for the pages actually registered in app.py.

Kept deliberately narrow: the template's original version walked a hard-coded
list of example pages, most of which quantms-web has deleted, so every one of
them failed with FileNotFoundError after a template sync.
"""

import ast
import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

# Pages AppTest.from_file can load in isolation. The rest call st.page_link,
# which needs the navigation context that only exists when app.py itself runs
# (loading them directly raises KeyError: 'url_pathname'), so they are covered
# indirectly by test_app_loads.
DIRECTLY_TESTABLE_PAGES = [
    "content/workflow_fileupload.py",
    "content/workflow_configure.py",
    "content/workflow_run.py",
]


def registered_pages() -> list[str]:
    """Every st.Page(Path("content", "...")) target declared in app.py."""
    tree = ast.parse(Path("app.py").read_text(encoding="utf-8"))
    pages = []
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "Page"
            and node.args
        ):
            target = node.args[0]
            if (
                isinstance(target, ast.Call)
                and isinstance(target.func, ast.Name)
                and target.func.id == "Path"
            ):
                parts = [a.value for a in target.args if isinstance(a, ast.Constant)]
                if parts:
                    pages.append("/".join(parts))
    return pages


def _init(apptest):
    with open("settings.json", "r", encoding="utf-8") as f:
        apptest.session_state.settings = json.load(f)
    apptest.session_state.settings["test"] = True
    apptest.secrets["workspace"] = "test"
    return apptest


@pytest.fixture
def launch(request):
    return _init(AppTest.from_file(request.param))


def test_registered_pages_exist():
    """Guard against a template sync registering pages this repo has deleted."""
    missing = [p for p in registered_pages() if not Path(p).is_file()]
    assert not missing, f"app.py registers pages that do not exist: {missing}"


@pytest.mark.parametrize("launch", DIRECTLY_TESTABLE_PAGES, indirect=True)
def test_page_loads(launch):
    launch.run(timeout=60)
    assert not launch.exception


def test_app_loads():
    app = _init(AppTest.from_file("app.py"))
    app.run(timeout=60)
    assert not app.exception
