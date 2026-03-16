import ast
from pathlib import Path
from streamlit.testing.v1 import AppTest
import pytest
import json


def get_pages_from_app():
    """Parse app.py AST to extract page paths from st.Page(Path(...)) calls."""
    tree = ast.parse(Path("app.py").read_text())
    pages = []
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "Page"
            and node.args
            and isinstance(node.args[0], ast.Call)
            and isinstance(node.args[0].func, ast.Name)
            and node.args[0].func.id == "Path"
        ):
            parts = [
                arg.value
                for arg in node.args[0].args
                if isinstance(arg, ast.Constant) and isinstance(arg.value, str)
            ]
            if parts:
                pages.append(str(Path(*parts)))
    return pages


def _uses_page_link(path: str) -> bool:
    """Return True if the file calls st.page_link(), which is incompatible with AppTest."""
    return "st.page_link(" in Path(path).read_text()


# Collect all content pages: those registered in app.py plus any other .py files
# in content/ (utility pages like digest.py, fragmentation.py, etc.).
# Exclude pages using st.page_link() — these require full st.navigation()
# context and cannot be launched in isolation via AppTest.
_app_pages = get_pages_from_app()
_all_content = sorted(
    str(p) for p in Path("content").glob("*.py") if p.name != "__init__.py"
)
_pages_to_test = sorted(
    p for p in set(_app_pages) | set(_all_content) if not _uses_page_link(p)
)


@pytest.fixture
def launch(request):
    test = AppTest.from_file(request.param)

    ## Initialize session state ##
    with open("settings.json", "r") as f:
        test.session_state.settings = json.load(f)
    test.session_state.settings["test"] = True
    test.secrets["workspace"] = "test"
    return test


# Test launching of all pages
@pytest.mark.parametrize("launch", _pages_to_test, indirect=True)
def test_launch(launch):
    """Test if all pages can be launched without errors."""
    launch.run(timeout=30)
    assert not launch.exception
