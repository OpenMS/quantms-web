"""Download Results Page."""
import streamlit as st
from src.common.common import page_setup
from src.common.guide import page_intro
from src.WorkflowTest import WorkflowTest
from src.common.combined_exports import build_missing, combined_paths, filtered_idxmls, load_run_params

params = page_setup()

wf = WorkflowTest()

st.title("Runs and exports" if st.session_state.location == "local" else "Runs and downloads")

page_intro(
    """
Every analysis you start is kept as its own run, so you can change a setting, run again
and compare. The list below shows all runs: 👁️ shows a run on the results pages (the
**Showing run** box at the top of each results page does the same), ✏️ renames it, 📦 packs its results,
settings and logs into one ZIP to download, and 🗑️ deletes it.

Below the list are the files of one run, singly or as one ZIP. The ones most people need:

- **quant_results/**: the protein and peptide abundance table (CSV, opens in Excel or R).
- **exports/**: all peptide-spectrum matches per run (`*_psms.tsv`) and each match's
  fragment ions (`*_fragment_ions.tsv`).
- **combined/**: the whole experiment in two files: every file's peptide-spectrum matches
  after the FDR filter in one idXML (`all_psms_filtered.idXML`, each match keeps its
  source file), and the experimental design as `design.tsv` (OpenMS format, which
  ProteomicsLFQ, MSstats and quantms read).
- **psm_filter/**, **percolator_results/**, **comet_results/**: identifications after each
  step in the open idXML format, one file per input file.
- The spectral library, if you built one.

The group comparison (fold changes and p-values) is downloaded from the table on the
**Compare groups** page with the table's download button.
""",
)



def combined_files_section() -> None:
    """The run's experiment-wide files, with a button to create them for older runs."""
    paths = combined_paths(wf.results_dir)
    missing = [k for k, p in paths.items() if not p.is_file()]
    if "idxml" in missing and not filtered_idxmls(wf.results_dir):
        return  # the run has no filtered identifications to merge
    with st.container(border=True):
        st.markdown("**Whole-experiment files**")
        st.markdown(
            f"- `combined/{paths['idxml'].name}`: all files' filtered peptide-spectrum matches in one idXML\n"
            f"- `combined/{paths['design'].name}`: the experimental design (file, condition, replicate, fraction)"
        )
        if not missing:
            st.caption("Both are in the list below, in the **combined** folder.")
            return
        st.caption("This run is older than these files. Create them from its results and settings:")
        if st.button("Create whole-experiment files", key="build-combined-exports"):
            with st.spinner("Merging the filtered PSMs with IDMerger..."):
                params = load_run_params(wf.runs.params_file(wf.run_id)) if wf.run_id else wf.params
                done = build_missing(wf.results_dir, params, wf.workflow_dir)
            if done.get("idxml") is False:
                st.error("IDMerger could not merge this run's idXML files.")
            else:
                st.rerun()


EXCLUDE = ["insight_cache"]
if wf.multi_run:
    st.markdown("#### Runs")
    wf.show_run_manager_section(exclude=EXCLUDE)
    st.markdown("#### Files of one run")
    if wf.select_run():
        combined_files_section()
        wf.ui.results_download_section(exclude=EXCLUDE, results_dir=wf.results_dir)
else:
    wf.show_results_download_section(exclude=EXCLUDE)
