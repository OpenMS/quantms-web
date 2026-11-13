"""How a run came about and the settings that differ from the defaults (src/viewer/history.py)."""
import json
from pathlib import Path

from src.viewer import history
from src.workflow.RunManager import RunManager

INI = """<?xml version="1.0" encoding="ISO-8859-1"?>
<PARAMETERS version="1.7.0">
  <NODE name="ProSE"><NODE name="1">
    <ITEM name="threads" value="1" type="int"/>
    <NODE name="Search">
      <ITEM name="enzyme" value="Trypsin" type="string"/>
      <ITEM name="min_size" value="7" type="int"/>
    </NODE>
    <ITEM name="score" value="NaN" type="double"/>
  </NODE></NODE>
</PARAMETERS>
"""


def _run(wf: Path, run_id: str, name: str, files, params, adds_to=None):
    folder = wf / "runs" / run_id
    (folder / "logs").mkdir(parents=True)
    run = {"name": name, "created": f"2026-11-1{run_id[-1]}T10:00:00", "status": "finished",
           "finished": "2026-11-13T10:05:00"}
    if adds_to:
        run["adds_to"] = adds_to
    (folder / "run.json").write_text(json.dumps(run))
    (folder / "params.json").write_text(json.dumps({"mzML-files": files}))
    step = {"tool": "ProSE", "files": [{"status": "done"}],
            "provenance": {"kind": "topp", "program": "ProSE", "params": params, "io": []}}
    (folder / "logs" / "progress.json").write_text(json.dumps({"steps": [step]}))


def test_timeline_and_changed_settings(tmp_path):
    wf = tmp_path / "topp-workflow"
    (wf / "ini").mkdir(parents=True)
    (wf / "ini" / "ProSE.ini").write_text(INI)
    _run(wf, "run-1", "Run 1", ["/data/a.mzML", "/data/b.mzML"],
         {"Search:enzyme": "Trypsin/P", "Search:min_size": "7", "threads": "4", "score": "nan"})
    _run(wf, "run-2", "More files", ["/data/c.mzML"], {"Search:enzyme": "Trypsin"}, adds_to="run-1")
    runs = RunManager(wf)

    events = history.timeline(runs, "run-1")
    assert [e["what"] for e in events] == [history.STARTED, history.FILES_ADDED]
    assert events[0]["files"] == ["a.mzML", "b.mzML"] and events[1]["name"] == "More files"
    assert "`c.mzML`" in history.event_markdown(events[1])

    # Defaults (also NaN in another spelling) and threads are left out; a value
    # that differed between the runs of the experiment lists both
    assert history.changed_settings(runs, "run-1") == [
        ("ProSE", "Search:enzyme", "Trypsin", "Trypsin/P | Trypsin"),
    ]


def test_run_without_report(tmp_path):
    wf = tmp_path / "topp-workflow"
    folder = wf / "runs" / "run-1"
    folder.mkdir(parents=True)
    (folder / "run.json").write_text(json.dumps({"name": "Imported", "status": "finished"}))
    runs = RunManager(wf)
    assert history.changed_settings(runs, "run-1") == []
    assert history.timeline(runs, "run-1")[0]["name"] == "Imported"
