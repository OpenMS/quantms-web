"""results_download as a Viewer component."""

import streamlit as st
from pathlib import Path
from src.WorkflowTest import WorkflowTest
from src.common.insight_tables import show_insight_table
from src.common.qpx_export import app_version, export_run, provenance_table, qpx_paths, read_provenance
from src.viewer.scope import scoped


def render(ctx) -> None:
    wf = WorkflowTest()
    wf.show_results_download_section(exclude=["insight_cache"])

    def qpx_section() -> None:
        """The run's QPX export: the provenance table, or a button to create the export for older runs."""
        results_dir = Path(wf.results_dir)
        if not (results_dir / "psm_filter").is_dir():
            return
        paths = qpx_paths(results_dir)
        st.markdown("#### Provenance (QPX)")
        if not all(p.is_file() for p in paths.values()):
            st.caption("This run is older than the QPX export.")
            if st.button("Create QPX files", key=scoped("build-qpx")):
                run_dir = results_dir.parent
                log = run_dir / "logs" / "commands-and-run-times.log"
                params = wf.runs.params(wf.run_id) if wf.run_id else wf.params
                with st.spinner("Writing the QPX files..."):
                    try:
                        export_run(
                            results_dir,
                            wf.workflow_dir,
                            params,
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
            name="qpx_provenance",
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

    if wf.multi_run and wf.run_id:
        qpx_section()
