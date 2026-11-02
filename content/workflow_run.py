import streamlit as st
from pathlib import Path
from src.common.common import page_setup
from src.common.guide import latest_empty_run, page_intro, show_empty_run
from src.common.run_estimate import estimate, format_minutes, total_memory
from src.WorkflowTest import WorkflowTest

params = page_setup()

wf = WorkflowTest()

st.title("Run the analysis")
page_intro(
    "Start the analysis and follow its progress. It runs on the server, so you can close "
    "this page and come back later with the same workspace. When it has finished, open "
    "**Run summary** under *Main results* in the sidebar.\n\n"
    "Each run is kept: to try other settings, change them and start again. The results "
    "pages show the newest run; switch runs with **Showing run** at the top of each results page."
)


@st.cache_data(show_spinner=False)
def _estimate(files: tuple, fasta: str, comet: tuple, threads: int, decoys: bool, mtimes: tuple, memory):
    """estimate(), cached on the inputs' modification times and the settings."""
    return estimate([Path(f) for f in files], Path(fasta) if fasta else None, dict(comet), threads, decoys, memory)


def show_estimate() -> None:
    """Time and memory the run would need, with advice when either is high."""
    run_params = wf.parameter_manager.get_parameters_from_json()
    selected = run_params.get("mzML-files") or []
    if isinstance(selected, str):
        selected = [selected]
    if not selected:
        return
    try:
        files = [str(f) for f in wf.file_manager.get_files(selected)]
        fasta = str(wf.file_manager.get_files([run_params["fasta-file"]])[0]) if run_params.get("fasta-file") else ""
    except (ValueError, KeyError):
        return
    instance = "CometAdapter" if run_params.get("analysis-mode", "LFQ") == "LFQ" else "CometAdapter-TMT"
    comet = wf.parameter_manager.get_merged_params(instance)
    comet = tuple(sorted((k, str(v)) for k, v in comet.items()))
    threads = wf.executor._get_max_threads()
    mtimes = tuple(Path(f).stat().st_mtime if Path(f).exists() else 0.0 for f in [*files, fasta] if f)
    est = _estimate(tuple(files), fasta, comet, threads, bool(run_params.get("generate-decoys", True)), mtimes, total_memory())

    with st.container(border=True):
        st.markdown("**Before you start: what this run needs**")
        c1, c2, c3 = st.columns(3)
        c1.metric("Spectra", f"{est.spectra:,}", help=f"In {est.files} file(s) (MS1 and MS2), searched against {est.proteins:,} proteins (before decoys).")
        c2.metric("Time", format_minutes(est.minutes), help=f"{est.parallel} file(s) are searched at the same time with {threads} thread(s).")
        memory = f"{est.peak_gb:.1f} GB" if est.peak_gb >= 1 else f"{est.peak_gb * 1024:.0f} MB"
        c3.metric(
            "Peak memory", memory,
            help="The most memory in use at one time: the files searched at once, or the quantification of the largest file.",
        )
        if est.total_gb:
            st.caption(
                f"{est.files} file(s), {est.parallel} searched at a time · this computer has {est.total_gb:.0f} GB of memory · "
                "a rough estimate from the file sizes, the database and the search settings."
            )
        if est.advice:
            text = "\n".join(f"- {a}" for a in est.advice)
            if est.level == "high":
                st.warning("**This run may need more memory than this computer has.**\n\n" + text, icon="⚠️")
            else:
                st.info("**This run will take a while.**\n\n" + text, icon="⏱️")
        with st.expander("How this is estimated"):
            st.markdown(
                "- **Memory**: Comet keeps the spectra of one batch in memory (about 0.3 MB per MS2 "
                "spectrum at high resolution; *spectrum_batch_size*, 20,000 by default, caps a batch), "
                "and the workflow searches as many files at once as it has threads. The quantification "
                "loads one mzML file at a time (about four times its size).\n"
                "- **Time**: Comet's speed measured on the demo data, scaled by the number of spectra, the "
                "database size (with decoys) and the search settings (semi-specific digestion, variable "
                "modifications, missed cleavages, precursor window), plus rescoring and quantification.\n\n"
                "The spectrum count comes from the mzML header (or the file size). Real runs can differ "
                "by a factor of two either way, more on a busy server."
            )


show_estimate()
# A run that kept no PSM or protein stops with an error; say why in plain words
if "workspace" in st.session_state:
    empty = latest_empty_run(st.session_state["workspace"])
    if empty is not None:
        _, run_name, diagnosis = empty
        show_empty_run(run_name, diagnosis, links=False)

wf.show_execution_section()

