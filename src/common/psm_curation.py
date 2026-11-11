"""Verdicts a user gives single PSMs: accepted or rejected, with a note.

Kept per run in ``results/curation.tsv``, so they travel with the run's ZIP
and outlive a session. A PSM is named by its file, scan, peptide and charge
rather than by a row number, so a verdict still finds its PSM when the run's
viewer caches are built again. Streamlit-free.
"""

import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import pandas as pd

CURATION_FILE = "curation.tsv"
VERDICTS = ("accepted", "rejected")
KEY = ["file", "scan", "sequence", "charge"]
COLUMNS = KEY + ["mz", "rt", "verdict", "note", "updated"]


def curation_path(results_dir: Path) -> Path:
    return Path(results_dir, CURATION_FILE)


def load_verdicts(results_dir: Path) -> pd.DataFrame:
    """Every verdict of the run, an empty table when there is none (or the file is unreadable)."""
    path = curation_path(results_dir)
    if not path.is_file():
        return pd.DataFrame(columns=COLUMNS)
    try:
        frame = pd.read_csv(path, sep="\t", dtype={"file": str, "sequence": str, "verdict": str, "note": str},
                            keep_default_na=False)
    except (OSError, ValueError, pd.errors.ParserError):
        return pd.DataFrame(columns=COLUMNS)
    for column in COLUMNS:
        if column not in frame:
            frame[column] = ""
    frame = frame[COLUMNS]
    frame["scan"] = pd.to_numeric(frame["scan"], errors="coerce").fillna(-1).astype(int)
    frame["charge"] = pd.to_numeric(frame["charge"], errors="coerce").fillna(0).astype(int)
    return frame[frame["verdict"].isin(VERDICTS) | (frame["note"] != "")].reset_index(drop=True)


def _match(frame: pd.DataFrame, psm: dict) -> pd.Series:
    return ((frame["file"] == str(psm["file"])) & (frame["scan"] == int(psm["scan"]))
            & (frame["sequence"] == str(psm["sequence"])) & (frame["charge"] == int(psm["charge"])))


def verdict_of(frame: pd.DataFrame, psm: dict) -> tuple[Optional[str], str]:
    """The verdict (or None) and the note of one PSM."""
    rows = frame[_match(frame, psm)]
    if rows.empty:
        return None, ""
    verdict = rows.iloc[-1]["verdict"]
    return (verdict if verdict in VERDICTS else None), str(rows.iloc[-1]["note"])


def set_verdict(results_dir: Path, psm: dict, verdict: Optional[str], note: str = "") -> pd.DataFrame:
    """Store the verdict and note of ``psm`` (file, scan, sequence, charge, mz, rt).

    A PSM left without a verdict and without a note is taken off the list.
    Returns the run's verdicts as they are now.
    """
    if verdict is not None and verdict not in VERDICTS:
        raise ValueError(f"verdict must be one of {VERDICTS} or None, not {verdict!r}")
    note = (note or "").replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()
    frame = load_verdicts(results_dir)
    frame = frame[~_match(frame, psm)]
    if verdict or note:
        row = {
            "file": str(psm["file"]), "scan": int(psm["scan"]), "sequence": str(psm["sequence"]),
            "charge": int(psm["charge"]), "mz": round(float(psm.get("mz") or 0.0), 6),
            "rt": round(float(psm.get("rt") or 0.0), 3), "verdict": verdict or "", "note": note,
            "updated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        }
        new = pd.DataFrame([row], columns=COLUMNS)
        frame = pd.concat([frame, new], ignore_index=True) if len(frame) else new
    _write(curation_path(results_dir), frame)
    return frame.reset_index(drop=True)


def _write(path: Path, frame: pd.DataFrame) -> None:
    """Write through a temporary file, so a reader never sees half a table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, tmp = tempfile.mkstemp(dir=path.parent, prefix=".curation-", suffix=".tsv")
    try:
        with os.fdopen(handle, "w", newline="") as out:
            frame.to_csv(out, sep="\t", index=False)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def verdicts_tsv(frame: pd.DataFrame) -> bytes:
    """The verdicts as a TSV download, accepted first."""
    order = {"accepted": 0, "rejected": 1}
    out = frame.assign(_o=frame["verdict"].map(order).fillna(2)).sort_values(["_o", "file", "scan"])
    return out.drop(columns="_o").to_csv(sep="\t", index=False).encode()
