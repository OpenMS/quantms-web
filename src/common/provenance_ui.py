"""A run's provenance from its QPX export: the table, or a button to create it for an older run."""
from pathlib import Path

import streamlit as st

from src.common.insight_tables import show_insight_table
from src.common.qpx_export import app_version, export_run, provenance_table, qpx_paths, read_provenance
from src.workflow.RunManager import RunManager


def has_provenance(runs: RunManager, run_id: str) -> bool:
    """True for a run that identified anything, so a QPX export exists or can be made."""
    return (runs.results_dir(run_id) / "psm_filter").is_dir()


def provenance(runs: RunManager, run_id: str) -> None:
    results_dir = Path(runs.results_dir(run_id))
    paths = qpx_paths(results_dir)
    if not all(p.is_file() for p in paths.values()):
        st.caption("This run is older than the QPX export.")
        if st.button("Create QPX files", key=f"build-qpx-{run_id}"):
            log = runs.run_dir(run_id) / "logs" / "commands-and-run-times.log"
            with st.spinner("Writing the QPX files..."):
                try:
                    export_run(
                        results_dir,
                        runs.workflow_dir,
                        runs.params(run_id),
                        log.read_text(encoding="utf-8", errors="ignore") if log.is_file() else "",
                        app_version(),
                    )
                except Exception as e:  # show the reason instead of a stack trace
                    st.error(f"Could not write the QPX files: {e}")
                    return
            st.rerun()
        return
    rows = read_provenance(paths["provenance"])
    table = provenance_table(rows)
    tools = [r for r in rows if r["tool_name"].startswith("OpenMS/")]
    st.caption(
        f"{len(tools)} OpenMS tools"
        + (f" (OpenMS {tools[0]['tool_version']})" if tools and tools[0].get("tool_version") else "")
        + f", {len(table)} parameters."
    )
    show_insight_table(
        table,
        name=f"qpx_provenance_{run_id}",
        cache_dir=results_dir / "insight_cache",
        source_files=[paths["provenance"]],
        title="provenance",
        height=420,
        layout="fitColumns",
        column_overrides={
            "step": {"width": 70},
            "tool": {"width": 160, "hozAlign": "left"},
            "version": {"width": 110, "hozAlign": "left"},
            "parameter": {"width": 180, "hozAlign": "left"},
            "value": {"hozAlign": "left"},
            "set by": {"width": 100, "hozAlign": "left"},
        },
    )
