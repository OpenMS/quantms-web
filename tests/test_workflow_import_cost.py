"""The workflow module is imported by every workflow page, so what it imports
at module level is paid on the first load of each of them. scipy.stats (about
0.6 s) is only needed by the group comparison and is imported there; mygene
(0.1 to 0.3 s on a cold disk) is used by nothing in the app."""

import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize("module", ["scipy.stats", "mygene"])
def test_import_leaves_module_unloaded(module):
    code = (
        "import sys\n"
        "import src.WorkflowTest\n"
        f"print({module!r} in sys.modules)\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", code],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip().splitlines()[-1] == "False"
