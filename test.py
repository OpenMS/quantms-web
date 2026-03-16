import unittest
import ast
import json
from pathlib import Path


def get_pages_from_app():
    """Parse app.py AST to extract page paths from st.Page(Path(...)) calls."""
    tree = ast.parse(Path("app.py").read_text(encoding="utf-8"))
    pages = []
    for node in ast.walk(tree):
        # Match st.Page(Path("content", "filename.py"), ...)
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


class TestSettingsJson(unittest.TestCase):
    def test_settings_json_exists(self):
        self.assertTrue(Path("settings.json").exists(), "settings.json file is missing")

    def test_settings_json_valid(self):
        with open("settings.json", "r") as f:
            settings = json.load(f)
        self.assertIn("app-name", settings)
        self.assertIn("version", settings)


class TestContentPagesExist(unittest.TestCase):
    def test_all_content_pages_exist(self):
        """Test that all content pages referenced by app.py exist."""
        pages = get_pages_from_app()
        self.assertTrue(
            len(pages) > 0,
            "No pages found in app.py — ensure app.py contains st.Page(Path(...)) calls",
        )
        for page in pages:
            self.assertTrue(Path(page).exists(), f"Content page {page} is missing")


if __name__ == '__main__':
    unittest.main()
