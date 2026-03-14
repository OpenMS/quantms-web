import unittest
import json
from pathlib import Path


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
        expected_pages = [
            "content/quickstart.py",
            "content/workflow_fileupload.py",
            "content/workflow_configure.py",
            "content/workflow_run.py",
            "content/results_database_search.py",
            "content/results_rescoring.py",
            "content/results_filtered.py",
            "content/results_abundance.py",
            "content/results_volcano.py",
            "content/results_pca.py",
            "content/results_heatmap.py",
            "content/results_library.py",
        ]
        for page in expected_pages:
            self.assertTrue(Path(page).exists(), f"Content page {page} is missing")


if __name__ == '__main__':
    unittest.main()
