"""Download Results Page."""
import streamlit as st
from pathlib import Path
from src.common.common import page_setup
from src.common.guide import page_intro
from src.WorkflowTest import WorkflowTest
from src.common.insight_tables import show_insight_table
from src.common.qpx_export import app_version, export_run, provenance_table, qpx_paths, read_provenance

params = page_setup()

wf = WorkflowTest()

st.title("Runs and exports" if st.session_state.location == "local" else "Runs and downloads")

page_intro(
    """
Every analysis you start is kept as its own run, so you can change a setting, run again
and compare. The list below shows all runs: 👁️ shows a run on the results pages (the
**Showing run** box at the top of each results page does the same), ✏️ renames it, 📦 packs its results,
settings and logs into one ZIP to download, and 🗑️ deletes it. The ZIP also holds the identified
spectra (**spectra/**), so after **Open a run from a ZIP** its results can be browsed without the mzML files.

Below the list are the files of one run, singly or as one ZIP. The ones most people need:

- **quant_results/**: the protein and peptide abundance table (CSV, opens in Excel or R).
- **exports/**: all peptide-spectrum matches per run (`*_psms.tsv`) and each match's
  fragment ions (`*_fragment_ions.tsv`).
- **psm_filter/**, **percolator_results/**, **comet_results/**: identifications after each
  step in the open idXML format, for OpenMS and other proteomics software.
- **qpx/**: the peptide-spectrum matches in the QPX format (`run.psm.parquet`, quantms
  Parquet, read by the `qpx` Python package, DuckDB, pandas or R arrow) with the run's
  provenance (`run.provenance.parquet`): every tool, its version and every parameter. The
  provenance is also shown below the file list.
- The spectral library, if you built one.

The group comparison (fold changes and p-values) is downloaded from the table on the
**Compare groups** page with the table's download button.
""",
)

wf.show_results_download_section(exclude=["insight_cache"])


def qpx_section() -> None:
    """The run's QPX export: the provenance table, or a button to create the export for older runs."""
    results_dir = Path(wf.results_dir)
    if not (results_dir / "psm_filter").is_dir():
        return
    paths = qpx_paths(results_dir)
    st.markdown("#### Provenance (QPX)")
    if not all(p.is_file() for p in paths.values()):
        st.caption(
            "This run is older than the QPX export. Create it from the run's results and logs: "
            "the filtered PSMs as `qpx/run.psm.parquet` and the provenance as `qpx/run.provenance.parquet`."
        )
        if st.button("Create QPX files", key="build-qpx"):
            run_dir = results_dir.parent
            log = run_dir / "logs" / "commands-and-run-times.log"
            params = wf.runs.params(wf.run_id) if wf.run_id else wf.params
            with st.spinner("Writing the QPX files..."):
                try:
                    export_run(
                        results_dir, wf.workflow_dir, params,
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
        f"How this run was made: {len(tools)} OpenMS tools"
        + (f" (OpenMS {tools[0]['tool_version']})" if tools and tools[0].get("tool_version") else "")
        + f", {len(table)} parameters. **set by** tells a value the workflow passed from one the tool "
        "kept as its default. Paths are relative to the workflow folder. The same table is in "
        "`qpx/run.provenance.parquet`."
    )
    show_insight_table(
        table, name="qpx_provenance", cache_dir=results_dir / "insight_cache",
        source_files=[paths["provenance"]], title="provenance", height=420,
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
