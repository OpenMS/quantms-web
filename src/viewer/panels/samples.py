"""results_samples as a Viewer component."""

from pathlib import Path
import streamlit as st
from src.common.guide import page_link
from src.common.insight_tables import show_insight_table, show_series
from src.common.results_helpers import get_abundance_data, get_sample_group_map, get_run_dir, get_run_params
from src.common.samples import COUNT_COLUMNS, channel_samples, condition_summary, lfq_samples, short_sample_names
from openms_insight.core.settings import setting_value
from src.viewer.scope import scoped, shared_state


def render(ctx) -> None:
    run_dir = get_run_dir(ctx.workspace)
    results_dir = run_dir / "results"
    cache_dir = results_dir / "insight_cache"
    workflow_params = get_run_params(ctx.workspace)
    analysis_mode = workflow_params.get("analysis-mode", "LFQ")

    def _mtimes(paths) -> tuple:
        return tuple(Path(p).stat().st_mtime if Path(p).exists() else 0.0 for p in paths)

    @st.cache_data(show_spinner="Counting per sample...")
    def _lfq_samples(workflow_dir: str, files: tuple, groups: tuple, mtimes: tuple):
        """lfq_samples, cached on the group assignments and the results' modification times."""
        return lfq_samples(workflow_dir, list(files), dict(groups))

    if analysis_mode == "LFQ":
        filter_dir = results_dir / "psm_filter"
        if not filter_dir.is_dir() or not any(filter_dir.glob("*.idXML")):
            st.info("No results available yet. Please run the workflow first.")
            page_link("run")
            return
        mzml_files = [str(f) for f in workflow_params.get("mzML-files", [])]
        group_map = {
            key[len("mzML-group-") :]: value
            for key, value in workflow_params.items()
            if key.startswith("mzML-group-") and value
        }
        sources = sorted(filter_dir.glob("*.idXML")) + sorted((results_dir / "quant_results").glob("*.csv"))
        samples = _lfq_samples(str(run_dir), tuple(mzml_files), tuple(sorted(group_map.items())), _mtimes(sources))
    else:
        abundance = get_abundance_data(ctx.workspace)
        if abundance is None:
            st.info("No quantification result available yet. Please run the workflow first.")
            page_link("run")
            return
        pivot_df, _, group_map = abundance
        sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
        # The report's first columns are the protein and its description; the channels follow
        sample_columns = [c for c in pivot_df.columns[4:] if c != "PeptideSequence"]
        samples = channel_samples(pivot_df, sample_columns, sample_group_map)
        sources = sorted((results_dir / "quant_results").glob("*.csv"))
    if samples.empty:
        st.info("No samples in the last run.")
        return

    def _table_height(rows: int) -> int:
        """Room for every row up to a page of 12, plus the title, toolbar, header and pager."""
        return 160 + 34 * min(rows, 12)

    counts = [c for c in COUNT_COLUMNS if samples[c].notna().any()]
    group_key = "|".join(f"{s}={c}" for s, c in zip(samples["Sample"], samples["Condition"]))
    state_manager = shared_state()
    # The count shown is a setting in the plot's toolbar (⚙)
    count_setting = {"identifier": "samples_metric", "label": "Count shown", "kind": "choice",
                     "options": counts, "default": counts[0]}
    metric = setting_value(state_manager, count_setting)
    if metric not in counts:
        metric = counts[0]
    bars = samples[["Sample", "Condition", metric]].copy()
    bars["Sample"] = bars["Sample"].map(short_sample_names(bars["Sample"]))
    # The plot beside the table of the samples, both of one height; the conditions below
    # (at half width the pager takes two lines)
    row_height = _table_height(len(samples)) + 30
    plot_column, table_column = st.columns(2, gap="medium")
    with plot_column:
        show_series(
            bars, f"samples_{metric.lower().replace(' ', '_')}",
            height=max(260, row_height), x_column="Sample", y_column=metric, series_column="Condition",
            mode="bars", orientation="h", x_label="Sample", y_label=metric, y_range_mode="nonnegative",
            state_manager=state_manager, toolbar_settings=[count_setting],
            note="Samples without a group are shown as Unassigned."
            if (samples["Condition"] == "Unassigned").any() else None,
        )
    with table_column:
        show_insight_table(
            samples[["Sample", "Condition"] + counts],
            name="samples_per_sample",
            cache_dir=cache_dir,
            source_files=sources,
            extra_key=group_key,
            title="Per sample",
            height=row_height,
            layout="fitColumns",
        )
    summary = condition_summary(samples)
    # Means go as text, since the table stores floats as float32 and would show 36.8 as 36.79999923706055
    for col in summary.columns[2:]:
        summary[col] = summary[col].map(lambda v: "" if v != v else f"{v:g}")
    show_insight_table(
        summary,
        name="samples_per_condition",
        cache_dir=cache_dir,
        source_files=sources,
        extra_key=group_key,
        title="Per condition",
        height=_table_height(len(summary)),
        layout="fitColumns",
    )
