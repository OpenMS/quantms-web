"""Shared page furniture that guides users through the app.

Every page opens with ``page_intro`` (what the page answers and, folded away,
the background a non-specialist needs to read it) and ends with ``next_steps``
(where to go from here), so the app reads as one path from upload to biology.
"""
import streamlit as st

from src.common.differential import STEPS, describe_settings, resolve_settings, run_differential
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map

# Page paths, so links survive a page being renamed in one place only.
PAGES = {
    "home": ("content/quickstart.py", "Home", "🏠"),
    "upload": ("content/workflow_fileupload.py", "Upload data", "📁"),
    "configure": ("content/workflow_configure.py", "Describe experiment", "⚙️"),
    "run": ("content/workflow_run.py", "Run analysis", "🚀"),
    "overview": ("content/results_overview.py", "Run summary", "🧾"),
    "proteins": ("content/results_proteins.py", "Identified proteins", "🧩"),
    "abundance": ("content/results_abundance.py", "Protein abundances", "📋"),
    "compare": ("content/compare_groups.py", "Compare groups", "⚖️"),
    "volcano": ("content/results_volcano.py", "Volcano plot", "🌋"),
    "pca": ("content/results_pca.py", "Sample similarity (PCA)", "📊"),
    "heatmap": ("content/results_heatmap.py", "Heatmap", "🔥"),
    "clustered": ("content/results_heatmap_clustered.py", "Clustered heatmap", "🧬"),
    "enrichment": ("content/enrichment.py", "GO enrichment", "📉"),
    "samples": ("content/results_samples.py", "Per-sample counts", "🧪"),
    "download": ("content/results_download.py", "Download results", "⬇️"),
    "filtering": ("content/filtering.py", "1. Filter proteins", "🧹"),
    "imputation": ("content/imputation.py", "2. Missing values", "🩹"),
    "normalization": ("content/normalization.py", "3. Normalize", "⚖️"),
    "statistics": ("content/statistical.py", "4. Statistical test", "🔢"),
    "qc": ("content/results_qc.py", "Identification quality", "🩺"),
    "filtered": ("content/results_filtered.py", "Peptide-spectrum matches", "🎯"),
    "search": ("content/results_database_search.py", "Raw search hits (Comet)", "🔬"),
    "rescoring": ("content/results_rescoring.py", "Rescored hits (Percolator)", "📈"),
    "mass_shifts": ("content/results_mass_shifts.py", "Mass shifts", "📐"),
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


def needs_results(message: str = "No results yet. Upload your data, describe the experiment and run the analysis first.") -> None:
    """Explain that the page needs a finished run, link to the start, and stop."""
    st.info(message)
    cols = st.columns(3)
    with cols[0]:
        page_link("upload")
    with cols[1]:
        page_link("configure")
    with cols[2]:
        page_link("run")
    st.stop()


def ensure_group_comparison(params: dict, show_errors: bool = True):
    """Make sure the group comparison exists, computing it with saved settings if not.

    Steps the user applied by hand on the Downstream Analysis pages are kept;
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
            st.markdown("Computed step by step on the Fine-tune pages.")
        st.caption("Change any step on the Fine-tune pages; the plots follow.")
        cols = st.columns(4)
        for col, key in zip(cols, ["filtering", "imputation", "normalization", "statistics"]):
            with col:
                page_link(key)
