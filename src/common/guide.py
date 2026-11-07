"""Shared page furniture that guides users through the app.

Every page opens with ``page_intro`` (what the page answers and, folded away,
the background a non-specialist needs to read it) and ends with ``next_steps``
(where to go from here), so the app reads as one path from upload to biology.
"""
from pathlib import Path

import streamlit as st

from src.common.differential import STEPS, describe_settings, resolve_settings, run_differential
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map

# Page paths, so links survive a page being renamed in one place only.
PAGES = {
    "home": ("content/quickstart.py", "Home", "🏠"),
    "configure": ("content/workflow_configure.py", "Describe experiment", "⚙️"),
    "run": ("content/workflow_run.py", "Run analysis", "🚀"),
    # Every result is a component of the results page
    "viewer": ("content/viewer.py", "Explore results", "👁️"),
    "export": ("content/workflow_export.py", "Export results", "📤"),
}


def page_link(key: str, label: str | None = None) -> None:
    path, default_label, icon = PAGES[key]
    try:
        st.page_link(path, label=label or default_label, icon=icon)
    except KeyError:
        # A page loaded on its own (AppTest.from_file) has no navigation to link into
        st.markdown(f"{icon} {label or default_label}")


def page_intro(lead: str, background: str | None = None, background_title: str = "Background") -> None:
    """The page's purpose in plain words, with optional background folded away."""
    st.markdown(lead)
    if background:
        with st.expander(f"📖 {background_title}"):
            st.markdown(background)


def next_steps(*steps: tuple[str, str]) -> None:
    """Footer linking to the pages that usually come next, each with a reason."""
    st.markdown("---")
    st.markdown("**Where to go next**")
    for key, why in steps:
        cols = st.columns([2, 5])
        with cols[0]:
            page_link(key)
        cols[1].markdown(why)


def needs_results(message: str = "No results yet. Describe your experiment (your runs, database and groups) and run the analysis first.") -> None:
    """Explain that the page needs a finished run, link to the start, and stop."""
    st.info(message)
    cols = st.columns(2)
    with cols[0]:
        page_link("configure")
    with cols[1]:
        page_link("run")
    st.stop()


def ensure_group_comparison(params: dict, show_errors: bool = True):
    """Make sure the group comparison exists, computing it with saved settings if not.

    Steps the user applied by hand in Edit comparison are kept;
    only the missing ones run. Returns the statistics table, or None when there
    is no abundance data or the samples do not form two groups.
    """
    workspace = st.session_state.get("workspace")
    result = get_abundance_data(workspace) if workspace else None
    if result is None:
        return st.session_state.get("statistics_df")
    pivot_df, _, group_map = result
    id_col = get_id_column(workspace, pivot_df)
    sample_group_map = get_sample_group_map(workspace, pivot_df, group_map)
    # Groups edited on Describe experiment apply without a new run, so a
    # comparison made with other groups is out of date
    previous = st.session_state.get("comparison_group_map")
    if previous is not None and previous != sample_group_map:
        for key in [*STEPS, "normalized_transform", "comparison_settings"]:
            st.session_state.pop(key, None)
    st.session_state["comparison_group_map"] = sample_group_map
    if st.session_state.get("statistics_df") is not None:
        return st.session_state["statistics_df"]
    settings = resolve_settings(params)
    try:
        with st.spinner("Comparing groups with the recommended settings..."):
            out = run_differential(
                pivot_df, sample_group_map, id_col, settings,
                done={k: st.session_state.get(k) for k in STEPS},
            )
    except ValueError as err:
        if show_errors:
            st.warning(str(err))
        return None
    for key in [*STEPS, "normalized_transform"]:
        st.session_state[key] = out[key]
    st.session_state["comparison_settings"] = describe_settings({**settings, "postproc-stat-test": out["stat_test"]})
    return out["statistics_df"]


def show_comparison_settings() -> None:
    """Fold-away note of how the group comparison was computed."""
    lines = st.session_state.get("comparison_settings")
    with st.expander("⚙️ How the groups were compared"):
        if lines:
            st.markdown("\n".join(f"- {line}" for line in lines))
        else:
            st.markdown("Computed step by step in Edit comparison.")
        st.caption("Change any step with **Edit comparison** above the run's components; the plots follow.")


@st.cache_data(show_spinner="Checking the run...")
def _empty_run_diagnosis(run_dir: str, mzml_files: tuple, params_json: str, mtimes: tuple):
    """diagnose_empty_run, cached on the run's result files."""
    import json

    from src.common.empty_result import diagnose_empty_run
    from src.common.overview import count_ms2_spectra

    def ms2_total():
        counts = [count_ms2_spectra(f) for f in mzml_files]
        known = [c for c in counts if c is not None]
        return sum(known) if known else None

    return diagnose_empty_run(Path(run_dir), json.loads(params_json), ms2_total)


def latest_empty_run(workspace) -> tuple[str, str, dict] | None:
    """(run id, run name, diagnosis) when the newest run stopped without results, else None.

    A run that keeps no PSM (or no protein) stops with an error, so the results pages,
    which list finished runs only, would show an older run or nothing at all.
    """
    import json

    from src.common.results_helpers import get_run_dir, get_run_params, get_runs
    from src.workflow.RunManager import FAILED

    runs = get_runs(workspace)
    newest = runs.list_runs()
    if not newest or newest[0].get("status") != FAILED:
        return None
    run = newest[0]
    run_dir = get_run_dir(workspace, run["id"])
    files = sorted((run_dir / "results").glob("*/*.idXML"))
    if not files:
        return None
    params = get_run_params(workspace, run["id"])
    mzml_files = tuple(str(f) for f in params.get("mzML-files", []))
    mtimes = tuple((str(f), f.stat().st_mtime) for f in files)
    diagnosis = _empty_run_diagnosis(
        str(run_dir), mzml_files, json.dumps(params, sort_keys=True, default=str), mtimes
    )
    if diagnosis is None:
        return None
    diagnosis["files"] = len(mzml_files)
    return run["id"], run.get("name") or run["id"], diagnosis


def show_empty_run(name: str, diagnosis: dict, links: bool = True) -> None:
    """Say which step emptied the run, at which FDR, and the likely causes."""
    if diagnosis["what"] == "protein":
        st.warning(
            f"**{name} identified no proteins.** Peptide-spectrum matches passed the PSM FDR, "
            "but the protein FDR filter left no protein, so there is nothing to quantify."
        )
    else:
        st.warning(
            f"**{name} found no peptide-spectrum matches.** The {diagnosis['step']} left no "
            "target PSM, so there are no peptides, proteins or quantities to show."
        )
    ms2 = diagnosis["ms2_spectra"]
    st.markdown(
        f"**Settings used:** PSM FDR {diagnosis['psm_fdr']} · protein FDR "
        f"{diagnosis['protein_fdr']}"
        + (f" · {ms2:,} MS2 spectra in {diagnosis['files']} file(s)" if ms2 is not None else "")
    )
    st.markdown(
        "**Target PSMs after each step:** "
        + " → ".join(
            f"{c['step']} **{c['targets']:,}{' ' + c['unit'] if c.get('unit') else ''}**"
            for c in diagnosis["counts"]
        )
    )
    st.markdown("**Likely causes**, most likely first:")
    st.markdown("\n".join(f"{i}. {cause}" for i, cause in enumerate(diagnosis["causes"], 1)))
    if links:
        cols = st.columns(3)
        with cols[0]:
            page_link("configure", "Change the settings")
        with cols[1]:
            page_link("run", "Run again")


def explain_empty_run(workspace, shown_run: str | None) -> None:
    """On a results page: explain a newest run that ended without results.

    With no finished run to show instead, the explanation is the page and it stops.
    Otherwise it is folded above the older run the page goes on to show.
    """
    found = latest_empty_run(workspace)
    if found is None:
        return
    _, name, diagnosis = found
    if shown_run is None:
        show_empty_run(name, diagnosis)
        st.stop()
    what = "no proteins" if diagnosis["what"] == "protein" else "no peptide-spectrum matches"
    with st.expander(f"⚠️ The newest run, {name}, ended with {what}. Why?"):
        show_empty_run(name, diagnosis)
