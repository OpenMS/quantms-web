"""Exporting (local app) or downloading (online) the files of a run.

Used by the Export page, which lists every run, and by the Export button of
Explore results, which offers the run on screen. A run's files are everything
in its folder (results, parameters, logs) except the visualization caches;
online they come as the run's ZIP, which **Open a run from a ZIP** reads back.
"""
from io import BytesIO
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import streamlit as st

from src.common.common import export_files, export_folder_input
from src.workflow.RunManager import RunManager

EXPORT_PAGE = "content/workflow_export.py"
# Folders that are caches of the viewer, not results
EXCLUDE = ["insight_cache"]


def is_local() -> bool:
    """True in the local app, which exports to a folder instead of downloading."""
    return st.session_state.get("location") == "local"


def verb() -> str:
    return "Export" if is_local() else "Download"


def run_files(runs: RunManager, run_id: str) -> list[tuple[Path, Path]]:
    """(file, path in the run's folder) of every file the run's export holds."""
    run_dir = runs.run_dir(run_id)
    files = []
    for path in sorted(run_dir.rglob("*")):
        rel = path.relative_to(run_dir)
        if path.is_file() and not any(part in EXCLUDE for part in rel.parts):
            files.append((path, rel))
    return files


def format_size(n_bytes: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n_bytes < 1024 or unit == "GB":
            return f"{n_bytes:.0f} {unit}" if unit == "B" else f"{n_bytes:.1f} {unit}"
        n_bytes /= 1024
    return ""


def total_size(files: list[tuple[Path, Path]]) -> int:
    return sum(f.stat().st_size for f, _ in files)


def safe_name(name: str) -> str:
    """A run name usable as a file or folder name."""
    text = str(name).replace(": ", " - ")
    cleaned = "".join("_" if c in '/\\:*?"<>|' or ord(c) < 32 else c for c in text).strip(" .")
    return cleaned or "run"


def _folders(runs: RunManager, run_ids: list[str]) -> dict[str, str]:
    """Run id to a folder name (its name, made unique)."""
    folders, used = {}, set()
    for run_id in run_ids:
        base = folder = safe_name((runs.get(run_id) or {}).get("name", run_id))
        n = 2
        while folder in used:
            folder, n = f"{base} ({n})", n + 1
        used.add(folder)
        folders[run_id] = folder
    return folders


def zip_runs(runs: RunManager, run_ids: list[str]) -> bytes:
    """One ZIP with a folder per run."""
    buffer = BytesIO()
    with ZipFile(buffer, "w", ZIP_DEFLATED) as zf:
        for run_id, folder in _folders(runs, run_ids).items():
            for path, rel in run_files(runs, run_id):
                zf.write(path, arcname=f"{folder}/{rel.as_posix()}")
    return buffer.getvalue()


def run_button(runs: RunManager, run_id: str, key: str, export_dir: Path | None = None,
               label: str | None = None, message_area=None) -> None:
    """One click downloads the run's ZIP (online) or copies its files to
    ``export_dir/<name>`` (local app), reporting the export in `message_area`."""
    name = (runs.get(run_id) or {}).get("name", run_id)
    if is_local():
        if st.button(label or "Export", key=f"export_{key}", icon=":material/drive_file_move:",
                     type="primary", width="stretch"):
            with message_area if message_area is not None else st.container():
                export_files(run_files(runs, run_id), Path(export_dir, safe_name(name)))
        return
    st.download_button(
        label or "Download ZIP",
        data=lambda: runs.zip(run_id, exclude=EXCLUDE).read_bytes(),
        file_name=f"{safe_name(name)}.zip", mime="application/zip",
        key=f"download_{key}", icon=":material/download:", type="primary",
        on_click="ignore", width="stretch",
    )


def all_runs_button(runs: RunManager, run_ids: list[str], export_dir: Path | None = None) -> None:
    """Export (or download) every run at once, a folder per run."""
    label = f"{verb()} all {len(run_ids)} runs" if len(run_ids) > 1 else f"{verb()} the run"
    if is_local():
        if st.button(label, key="export_all_runs", icon=":material/drive_file_move:", type="primary"):
            for run_id, folder in _folders(runs, run_ids).items():
                export_files(run_files(runs, run_id), Path(export_dir, folder))
        return
    app = safe_name(st.session_state.get("settings", {}).get("app-name", "results"))
    st.download_button(
        label, data=lambda: zip_runs(runs, run_ids), file_name=f"{app}-results.zip",
        mime="application/zip", key="download_all_runs", icon=":material/download:",
        type="primary", on_click="ignore",
    )


def file_table(files: list[tuple[Path, Path]]) -> None:
    """The files an export holds, with their sizes."""
    st.dataframe(
        {"File": [rel.as_posix() for _, rel in files],
         "Size": [format_size(f.stat().st_size) for f, _ in files]},
        hide_index=True, width="stretch",
    )


def _dialog_body(runs: RunManager, run_id: str) -> None:
    name = (runs.get(run_id) or {}).get("name", run_id)
    files = run_files(runs, run_id)
    st.markdown(f"**{name}**: {len(files)} files, {format_size(total_size(files))}")
    file_table(files)
    export_dir = export_folder_input("viewer-export") if is_local() else None
    run_button(runs, run_id, "viewer_dialog", export_dir)
    try:
        st.page_link(EXPORT_PAGE, label=f"All runs on the {verb()} page", icon=":material/folder_zip:")
    except KeyError:
        # A page loaded on its own (AppTest.from_file) has no navigation to link into
        pass


@st.dialog("Download results", width="large")
def _download_dialog(runs: RunManager, run_id: str) -> None:
    _dialog_body(runs, run_id)


@st.dialog("Export results", width="large")
def _export_dialog(runs: RunManager, run_id: str) -> None:
    _dialog_body(runs, run_id)


def open_dialog(runs: RunManager, run_id: str | None) -> None:
    """The shown run's export (Explore results), or the Export page without one."""
    if run_id is None:
        st.switch_page(EXPORT_PAGE)
    (_export_dialog if is_local() else _download_dialog)(runs, run_id)


def viewer_button(runs: RunManager, run_id: str | None) -> None:
    """The experiment row's button on Explore results: the shown run's export."""
    if st.button(verb(), icon=":material/download:" if not is_local() else ":material/output:",
                 key="viewer_export", type="primary", width="content",
                 help=f"{verb()} the files of the run shown below"):
        open_dialog(runs, run_id)
