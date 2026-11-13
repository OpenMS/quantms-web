"""How a run came about and the settings it was made with (Streamlit-free).

``timeline`` lists one event per start of a run's experiment (the run that
started it, then every run that added files to it), ``changed_settings`` the
TOPP tool settings those starts used that differ from the tools' defaults.
"""
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from src.workflow.processing_report import report_for
from src.workflow.parameter_names import display_name
from src.workflow.RunManager import FOLLOW_UP_KEY, RunManager

FILES_KEY = "mzML-files"
STARTED, FILES_ADDED = "Started", "Files added"
# Parameters that name files or threads, not settings
NOT_SETTINGS = {"threads", "in", "out", "database", "fasta", "design", "out_mzTab", "out_msstats",
                "out_cxml", "out_triqler", "in_features", "id", "out_id"}


def _time(value: Optional[str]) -> str:
    try:
        return datetime.fromisoformat(value).strftime("%Y-%m-%d %H:%M")
    except (TypeError, ValueError):
        return value or ""


def _files(runs: RunManager, run_id: str) -> List[str]:
    value = runs.params(run_id).get(FILES_KEY) or []
    return [Path(f).name for f in ([value] if isinstance(value, str) else value)]


def timeline(runs: RunManager, run_id: str) -> List[dict]:
    """The starts of ``run_id``'s experiment, oldest first: ``what``, ``name``,
    ``when`` (text), ``files`` (names) and ``detail``."""
    events = []
    for index, member in enumerate(runs.experiment(run_id) or [run_id]):
        run = runs.get(member) or {}
        follow = runs.params(member).get(FOLLOW_UP_KEY)
        detail = ""
        if index == 0 and isinstance(follow, dict) and follow.get("name"):
            detail = f"with the settings of {follow['name']}"
        events.append({
            "what": STARTED if index == 0 else FILES_ADDED,
            "name": run.get("name", member),
            "when": _time(run.get("created")),
            "files": _files(runs, member),
            "detail": detail,
        })
    return events


def event_markdown(event: dict) -> str:
    """One event as a line of Markdown."""
    text = f"**{event['what']}** {event['when']} as *{event['name']}*".replace("  ", " ")
    names = event["files"]
    if names:
        shown = ", ".join(f"`{n}`" for n in names[:4])
        if len(names) > 4:
            shown += f" and {len(names) - 4} more"
        text += f", on {shown}"
    if event["detail"]:
        text += f" ({event['detail']})"
    return text


def ini_defaults(ini_path: Path) -> Dict[str, str]:
    """parameter path -> default value (text) of a TOPP tool's .ini, ``{}`` without one."""
    try:
        root = ET.parse(ini_path).getroot()
    except (ET.ParseError, OSError):
        return {}
    found: Dict[str, str] = {}

    def tag(el) -> str:
        return el.tag.rsplit("}", 1)[-1] if isinstance(el.tag, str) else ""

    def walk(el, parts):
        for child in el:
            kind = tag(child)
            if kind == "NODE":
                walk(child, parts + [child.get("name") or ""])
            elif kind == "ITEM":
                found[":".join(parts[2:] + [child.get("name") or ""])] = child.get("value") or ""
            elif kind == "ITEMLIST":
                values = [i.get("value") or "" for i in child if tag(i) == "LISTITEM"]
                found[":".join(parts[2:] + [child.get("name") or ""])] = ", ".join(values)

    # <PARAMETERS><NODE name="Tool"><NODE name="1"> ...: the path starts below those two
    walk(root, [])
    return found


def _same(value: str, default: str) -> bool:
    if value == default or value.strip().lower() == default.strip().lower():
        return True
    try:
        return float(value) == float(default)
    except (TypeError, ValueError):
        return False


def changed_settings(runs: RunManager, run_id: str) -> List[Tuple[str, str, str, str]]:
    """(tool, parameter, default, value) of the TOPP tool settings the runs of
    ``run_id``'s experiment used that differ from the tools' defaults; a value
    that differed between the runs lists each one."""
    per_tool: Dict[str, Dict[str, List[str]]] = {}
    for member in runs.experiment(run_id) or [run_id]:
        try:
            steps = report_for(runs.workflow_dir, member, settings={})["steps"]
        except Exception:  # an imported run with damaged logs has no report
            steps = []
        for step in steps:
            if step.get("kind") != "topp":
                continue
            values = per_tool.setdefault(step.get("program") or step.get("tool"), {})
            for key, value in (step.get("params") or {}).items():
                found = values.setdefault(key, [])
                if str(value) not in found:
                    found.append(str(value))
    rows = []
    for tool, values in per_tool.items():
        defaults = ini_defaults(Path(runs.workflow_dir, "ini", f"{tool}.ini"))
        for key, found in values.items():
            if key.split(":")[-1] in NOT_SETTINGS:
                continue
            default = defaults.get(key)
            if default is not None and len(found) == 1 and _same(found[0], default):
                continue
            entry = display_name(tool, key)
            label = entry["name"] if entry else key
            rows.append((tool, label, "" if default is None else default, " | ".join(found)))
    return rows
