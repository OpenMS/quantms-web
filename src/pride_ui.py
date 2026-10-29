"""File Upload tab that imports a PRIDE Archive dataset for reprocessing.

The user looks up an accession, picks the MS runs (mzML, or Thermo RAW which
FileConverter turns into mzML), FASTA and SDRF files, and downloads them in a
background process into the workflow's input folders. An SDRF file then sets
the analysis mode, sample groups and search settings on the Configure page.
"""
import json
import multiprocessing
import time
from pathlib import Path

import pandas as pd
import streamlit as st

from src import pride
from src.common.insight_tables import show_insight_table

LOOKUP_KEY = "pride-lookup"
STARTED_KEY = "pride-job-started"
# Downloads are capped on public servers, where the workspace volume is shared
DEFAULT_ONLINE_LIMIT_GB = 50


def _dirs(wf) -> tuple[Path, Path, Path]:
    inputs = Path(wf.workflow_dir, "input-files")
    return Path(wf.workflow_dir, "pride-import"), Path(inputs, "mzML-files"), Path(inputs, "fasta-file")


def _download_limit_bytes() -> int | None:
    if st.session_state.get("location") != "online":
        return None
    limit = st.session_state.get("settings", {}).get("pride_import", {}).get(
        "max_download_gb_online", DEFAULT_ONLINE_LIMIT_GB
    )
    return int(float(limit) * 1e9) if limit else None


def _format_size(n_bytes: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n_bytes < 1000:
            return f"{n_bytes:.1f} {unit}" if unit != "B" else f"{int(n_bytes)} B"
        n_bytes /= 1000
    return f"{n_bytes:.1f} TB"


def _lookup(accession_text: str) -> None:
    try:
        accession = pride.normalize_accession(accession_text)
    except ValueError as e:
        st.error(str(e))
        return
    try:
        with st.spinner(f"Fetching {accession} from PRIDE..."):
            project = pride.fetch_project(accession)
            files = pride.fetch_files(accession)
    except LookupError as e:
        st.error(str(e))
        return
    except Exception as e:
        st.error(f"Could not reach the PRIDE Archive: {e}")
        return
    st.session_state[LOOKUP_KEY] = {"project": project, "files": pride.as_dicts(files)}


def _show_project(project: dict) -> None:
    st.markdown(f"#### [{project['accession']}]({project['url']}): {project['title']}")
    facts = {
        "Organisms": project["organisms"],
        "Instruments": project["instruments"],
        "Modifications": project["ptms"],
        "Quantification": project["quantification"],
        "Experiment type": project["experiment_types"],
    }
    cols = st.columns(len(facts))
    for col, (label, values) in zip(cols, facts.items()):
        col.markdown(f"**{label}**  \n{', '.join(values) if values else 'not given'}")
    with st.expander("Description and protocols"):
        for label, key in (("Description", "description"),
                           ("Sample processing protocol", "sample_protocol"),
                           ("Data processing protocol", "data_protocol")):
            if project[key]:
                st.markdown(f"**{label}**")
                st.write(project[key])
        for reference in project["references"]:
            st.caption(reference)


def _run_stem(name: str) -> str:
    return name.lower().removesuffix(".gz").rsplit(".", 1)[0]


def _default_ms_files(files: list[dict]) -> list[str]:
    """Every run once: its mzML if the dataset has one, else its RAW file."""
    mzml = [f["name"] for f in files if f["kind"] == pride.KIND_MZML]
    stems = {_run_stem(n) for n in mzml}
    raw = [f["name"] for f in files if f["kind"] == pride.KIND_RAW and _run_stem(f["name"]) not in stems]
    return mzml + raw


def _start_import(wf, accession: str, selected: list[dict], keep_raw: bool, fallback: dict) -> None:
    job_dir, mzml_dir, fasta_dir = _dirs(wf)
    job_dir.mkdir(parents=True, exist_ok=True)
    example_files = {
        str(mzml_dir): [Path(f).name for f in fallback.get("mzML-files", [])],
        str(fasta_dir): [Path(f).name for f in fallback.get("fasta-file", [])],
    }
    pride.status_path(job_dir).unlink(missing_ok=True)
    Path(job_dir, "import.log").unlink(missing_ok=True)
    process = multiprocessing.Process(
        target=pride.run_import,
        args=(str(job_dir), accession, selected, str(mzml_dir), str(fasta_dir), keep_raw,
              example_files),
        daemon=False,
    )
    process.start()
    st.session_state[STARTED_KEY] = time.time()


def _file_selection(wf, lookup: dict, fallback: dict) -> None:
    project, files = lookup["project"], lookup["files"]
    by_name = {f["name"]: f for f in files}

    st.markdown("##### Dataset files")
    table = pd.DataFrame(files)[["name", "kind", "category", "size"]]
    table["size"] = table["size"].map(_format_size)
    table = table.rename(columns={"name": "File", "kind": "Type", "category": "PRIDE category",
                                  "size": "Size"})
    show_insight_table(table, name="pride-files", extra_key=project["accession"],
                       title=f"{project['accession']}_files", height=300, page_size=25)

    ms_options = [f["name"] for f in files if f["kind"] in (pride.KIND_MZML, pride.KIND_RAW)]
    fasta_options = [f["name"] for f in files if f["kind"] == pride.KIND_FASTA]
    sdrf_options = [f["name"] for f in files if f["kind"] == pride.KIND_SDRF]
    if not ms_options:
        st.warning("This dataset has no mzML or Thermo RAW files that can be imported.")
        return

    # Widget keys carry the accession so a new lookup starts from its defaults
    acc = project["accession"]
    ms_files = st.multiselect(
        "MS runs to import", ms_options, default=_default_ms_files(files), key=f"pride-ms-{acc}",
        help="mzML files are used as they are; Thermo RAW files are converted to mzML with FileConverter.",
    )
    c1, c2 = st.columns(2)
    fasta_files = c1.multiselect(
        "FASTA databases", fasta_options, default=fasta_options, key=f"pride-fasta-{acc}",
        help="Most submissions do not include their database. Without one, add a FASTA in the "
             "'FASTA database' tab.",
        disabled=not fasta_options,
        placeholder="None in this dataset" if not fasta_options else "Choose files",
    )
    sdrf_files = c2.multiselect(
        "SDRF sample metadata", sdrf_options, default=sdrf_options, key=f"pride-sdrf-{acc}",
        help="An SDRF file describes each run's samples and search settings; it can set up the "
             "workflow after the import.",
        disabled=not sdrf_options,
        placeholder="None in this dataset" if not sdrf_options else "Choose files",
    )

    selected = [by_name[n] for n in ms_files + fasta_files + sdrf_files]
    has_raw = any(by_name[n]["kind"] == pride.KIND_RAW for n in ms_files)
    keep_raw = False
    can_convert, reason = pride.raw_conversion_available()
    if has_raw:
        if can_convert:
            keep_raw = st.checkbox("Keep RAW files after conversion", value=False,
                                   key=f"pride-keep-raw-{acc}")
        else:
            st.error(f"Thermo RAW files cannot be converted here: {reason} "
                     "Choose mzML files, or run the app from the Docker image.")

    total = sum(f["size"] for f in selected)
    limit = _download_limit_bytes()
    too_large = limit is not None and total > limit
    st.caption(f"{len(selected)} file(s), {_format_size(total)} to download.")
    if too_large:
        st.error(f"This server imports at most {_format_size(limit)} per dataset. "
                 "Choose fewer runs.")

    if st.button("Download and import", type="primary", use_container_width=True,
                 disabled=not ms_files or too_large or (has_raw and not can_convert),
                 key="pride-start"):
        _start_import(wf, acc, selected, keep_raw, fallback)
        st.rerun()


def _progress(job_dir: Path) -> None:
    status = pride.read_status(job_dir)
    if not status:
        st.info("Starting the import...")
        return
    done, total = status.get("bytes_done", 0), status.get("bytes_total", 0) or 1
    current = status.get("current", "")
    step = {"download": f"downloading {current}", "decompress": f"decompressing {current}",
            "convert": f"converting {current} to mzML"}.get(status.get("step"), current)
    number = min(status.get("done_files", 0) + 1, status.get("total_files", 0))
    st.progress(
        min(done / total, 1.0),
        text=f"{status.get('accession')}: file {number} of {status.get('total_files')}, {step}"
             f" ({_format_size(done)} of {_format_size(total)} downloaded)",
    )
    if st.button("Cancel import", key="pride-cancel"):
        pride.cancel(job_dir)
        time.sleep(1)
        st.rerun()


def _running(job_dir: Path) -> bool:
    multiprocessing.active_children()  # reap a finished import so its pid is gone
    if pride.is_running(job_dir):
        return True
    # The process may not have written its first status yet
    started = st.session_state.get(STARTED_KEY, 0)
    return not pride.read_status(job_dir) and time.time() - started < 30


@st.fragment(run_every=2)
def _progress_fragment(job_dir: Path) -> None:
    if not _running(job_dir):
        st.rerun(scope="app")
    _progress(job_dir)


def _show_result(job_dir: Path) -> dict:
    status = pride.read_status(job_dir)
    if not status:
        return status
    state = status.get("state")
    imported, errors = status.get("imported", []), status.get("errors", [])
    message = f"Import of {status.get('accession')}: {len(imported)} file(s) imported"
    if errors:
        message += f", {len(errors)} failed"
    if state == "cancelled":
        st.warning(f"{message}. The import was cancelled.")
    elif errors:
        st.warning(f"{message}.")
    else:
        st.success(f"{message}.")
    if errors:
        with st.expander("Errors"):
            for error in errors:
                st.write(error)
    log = Path(job_dir, "import.log")
    if log.exists():
        with st.expander("Import log"):
            st.code(log.read_text(encoding="utf-8", errors="replace")[-20000:], language=None)
    return status


def apply_to_workflow(wf, mzml_files: list[str], fasta_files: list[str], settings: dict | None) -> None:
    """Select the imported files and apply SDRF settings to the workflow's parameters.

    Same pattern as presets: write params.json, then drop the widgets' session
    keys so the Configure page re-initializes them from the file.
    """
    pm = wf.parameter_manager
    params = pm.get_parameters_from_json()
    drop = []

    def set_param(key, value):
        params[key] = value
        drop.append(f"{pm.param_prefix}{key}")

    def set_topp(tool, key, value):
        params.setdefault(tool, {})[key] = value
        session_key = f"{pm.topp_param_prefix}{tool}:1:{key}"
        drop.extend([session_key, f"{session_key}_display"])

    if mzml_files:
        set_param("mzML-files", mzml_files)
    if fasta_files:
        set_param("fasta-file", fasta_files[0])
    if settings:
        set_param("analysis-mode", settings["mode"])
        for name, group in settings["groups"].items():
            set_param(f"mzML-group-{name}", group)
        if settings["tmt_type"]:
            set_topp("IsobaricAnalyzer-TMT", "type", settings["tmt_type"])
        for index, group in settings["tmt_groups"].items():
            set_param(f"TMT-group-sample{int(index) + 1}", group)
        comet = "CometAdapter-TMT" if settings["mode"] == "TMT" else "CometAdapter"
        for key, value in settings["comet"].items():
            set_topp(comet, key, value)

    with open(pm.params_file, "w", encoding="utf-8") as f:
        json.dump(params, f, indent=4)
    for key in drop:
        st.session_state.pop(key, None)
    wf.params = params
    wf.ui.params = params


def _setup_section(wf, job_dir: Path, status: dict) -> None:
    _, mzml_dir, _ = _dirs(wf)
    imported = [Path(p) for p in status.get("imported", [])]
    mzml_files = [str(p) for p in imported if p.suffix.lower() == ".mzml" and p.exists()]
    fasta_files = [str(p) for p in imported
                   if pride.file_kind(p.name) == pride.KIND_FASTA and p.exists()]
    if not mzml_files:
        return

    st.markdown("##### Set up the workflow")
    sdrf_dir = Path(job_dir, "sdrf")
    uploaded = st.file_uploader(
        "Add your own SDRF file", type=["tsv", "txt"], key="pride-sdrf-upload",
        help="Use this when the dataset has no SDRF, or to replace it with a corrected one.",
    )
    if uploaded is not None:
        sdrf_dir.mkdir(parents=True, exist_ok=True)
        Path(sdrf_dir, uploaded.name).write_bytes(uploaded.getbuffer())

    sdrf_options = sorted(str(p) for p in sdrf_dir.glob("*")) if sdrf_dir.exists() else []
    sdrf_choice = st.selectbox(
        "Sample metadata (SDRF)", ["Don't use an SDRF file"] + sdrf_options,
        index=1 if sdrf_options else 0, key="pride-sdrf-choice",
        format_func=lambda p: Path(p).name if p in sdrf_options else p,
    )

    settings = None
    if sdrf_choice in sdrf_options:
        try:
            settings = pride.sdrf_settings(pride.read_sdrf(sdrf_choice),
                                           [Path(p).name for p in mzml_files])
        except (OSError, ValueError) as e:
            st.error(f"Could not read {Path(sdrf_choice).name}: {e}")
        if settings:
            preview = pd.DataFrame(pride.settings_preview(settings))
            show_insight_table(preview, name="pride-sdrf-settings", extra_key=sdrf_choice,
                               title="SDRF settings", height=300, page_size=25)
            for note in settings["notes"]:
                st.caption(note)

    if not fasta_files:
        st.warning("No FASTA database was imported. Add the protein database for this dataset's "
                   "organism in the 'FASTA database' tab.")
    summary = f"Selects the {len(mzml_files)} imported mzML file(s)"
    if fasta_files:
        summary += f" and {Path(fasta_files[0]).name}"
    if settings:
        summary += ", and applies the SDRF settings above"
    st.caption(summary + ". Review them on the Describe experiment page before running.")
    if st.button("Use for the workflow", type="primary", use_container_width=True,
                 key="pride-apply"):
        apply_to_workflow(wf, mzml_files, fasta_files, settings)
        st.success("Done. Open the Describe experiment page to review the settings.")


def pride_import_section(wf, fallback: dict = None) -> None:
    """The PRIDE import tab.

    ``fallback`` maps an upload key ("mzML-files", "fasta-file") to the example
    files its upload tab adds; they are cleared when the import replaces them.
    """
    st.markdown(
        "Reprocess a public dataset from the [PRIDE Archive](https://www.ebi.ac.uk/pride/archive): "
        "enter its accession to download the runs into this workspace."
    )
    job_dir, _, _ = _dirs(wf)

    if _running(job_dir):
        _progress_fragment(job_dir)
        return

    with st.form("pride-lookup-form", border=False):
        c1, c2 = st.columns([3, 1])
        accession = c1.text_input("PRIDE accession", placeholder="PXD000001",
                                  label_visibility="collapsed")
        looked_up = c2.form_submit_button("Look up", use_container_width=True)
    if looked_up:
        _lookup(accession)

    lookup = st.session_state.get(LOOKUP_KEY)
    if lookup:
        _show_project(lookup["project"])
        _file_selection(wf, lookup, fallback or {})

    status = _show_result(job_dir)
    if status and status.get("state") in ("done", "cancelled", "failed"):
        _setup_section(wf, job_dir, status)
