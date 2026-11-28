"""OpenDDA offers Workspace Sync to every user, for studies too large for the browser upload."""

import json
from pathlib import Path

SETTINGS = json.loads(Path(__file__).resolve().parents[1].joinpath("settings.json").read_text())


def test_workspace_sync_is_on_for_everyone():
    sync = SETTINGS["workspace_sync"]
    assert sync["enabled"] is True
    assert sync["developer_only"] is False


def test_workspace_sync_runs_the_workflow_class():
    """The cloud-run mapping names the workflow folder and class that exist."""
    sync = SETTINGS["workspace_sync"]
    (folder, target), = sync["workflows"].items()
    module, cls = target.split(":")
    assert folder == "topp-workflow"
    path = Path(__file__).resolve().parents[1] / (module.replace(".", "/") + ".py")
    assert path.exists() and f"class {cls}(" in path.read_text()


def test_app_registers_the_sync_page():
    """The sidebar panel links content/workspace_sync.py; without the page every page fails."""
    app = Path(__file__).resolve().parents[1].joinpath("app.py").read_text()
    assert 'Path("content", "workspace_sync.py")' in app
    assert Path(__file__).resolve().parents[1].joinpath("content", "workspace_sync.py").exists()
