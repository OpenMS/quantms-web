"""The app's wording of its tools in the methods section (settings.json ``methods``)."""

import json
import re
from pathlib import Path

from src.workflow.methods_text import methods_text
from src.workflow.processing_report import build_report, command_io
from src.workflow.progress import RunProgress, read_progress

SETTINGS = json.loads(Path("settings.json").read_text(encoding="utf-8"))


def _tools_run():
    """Every TOPP tool and instance name the workflow passes to run_topp()."""
    source = Path("src/WorkflowTest.py").read_text(encoding="utf-8")
    tools = set(re.findall(r'run_topp\(\s*"([^"]+)"', source))
    tools |= set(re.findall(r'tool_instance_name="([^"]+)"', source))
    return tools


def _covered(name, wording):
    # methods_text looks a step up by its instance name, then by its tool's name
    return name in wording or any(name.startswith(f"{tool}-") for tool in wording)


def test_every_tool_the_workflow_runs_is_worded():
    wording = SETTINGS["methods"]["tools"]
    missing = sorted(t for t in _tools_run() if not _covered(t, wording))
    assert not missing, f"no methods wording for {missing}"
    for entry in wording.values():
        assert entry["text"][0].isupper() and not entry["text"].endswith(".")


def _step(progress, tool, io, params=None, program=None):
    step = progress.start_step(
        tool,
        [Path(c["inputs"][0]).name for c in io],
        provenance={"kind": "topp", "program": program or tool, "version": "3.6.0",
                    "params": params or {}, "io": io},
    )
    for i in range(len(io)):
        progress.file_started(step, i)
        progress.file_finished(step, i, True, returncode=0)
    progress.end_step(step)


def test_a_label_free_run_reads_as_methods(tmp_path):
    progress = RunProgress(tmp_path)
    _step(progress, "DecoyDatabase", command_io({"in": [["/d/human.fasta"]], "out": ["/d/human.decoy.fasta"]}, 1))
    _step(progress, "ProSE", command_io(
        {"in": ["/d/A.mzML", "/d/B.mzML"], "database": ["/d/human.decoy.fasta"],
         "out": ["/r/A.idXML", "/r/B.idXML"]}, 2))
    _step(progress, "PercolatorAdapter", command_io(
        {"in": ["/r/A.idXML", "/r/B.idXML"], "out": ["/r/A_perc.idXML", "/r/B_perc.idXML"]}, 2))
    _step(progress, "IDFilter", command_io(
        {"in": ["/r/A_perc.idXML", "/r/B_perc.idXML"], "out": ["/r/A_f.idXML", "/r/B_f.idXML"]}, 2),
        params={"score:psm": "0.01"})
    _step(progress, "ProteomicsLFQ", command_io(
        {"in": [["/d/A.mzML", "/d/B.mzML"]], "ids": [["/r/A_f.idXML", "/r/B_f.idXML"]],
         "out": ["/r/out.mzTab"]}, 1))
    report = build_report(read_progress(tmp_path), app_name="OpenDDA", app_version="1.3.0", workflow="WF")

    text = methods_text(report, SETTINGS, names=lambda *a: None)

    assert "Then, spectra were searched against the protein database with ProSE on each file" in text
    assert "peptide-spectrum matches were rescored with Percolator [4] on each file" in text
    assert "peptides and proteins were quantified label-free with ProteomicsLFQ on all files together" in text
    assert "Elias JE, Gygi SP." in text and "Käll L," in text
