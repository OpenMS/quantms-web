import json


def test_settings_json_valid():
    """Test that settings.json exists and contains required fields."""
    with open("settings.json", "r") as f:
        settings = json.load(f)
    assert "app-name" in settings
    assert "version" in settings


def test_content_pages_exist():
    """Test that all content pages referenced by app.py exist."""
    from pathlib import Path

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
        assert Path(page).exists(), f"Content page {page} is missing"
