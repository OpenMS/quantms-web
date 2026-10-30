import os
# Polars' default (CPU-core-count-sized) native thread pool access-violation
# crashes the whole process on some high-core-count Windows machines when
# invoked from Streamlit's script-runner thread. Must be set before polars
# is imported anywhere (including transitively via openms_insight).
os.environ.setdefault("POLARS_MAX_THREADS", "1")

import streamlit as st
from pathlib import Path
import json
# For some reason the windows version only works if this is imported here
import pyopenms

if "settings" not in st.session_state:
        with open("settings.json", "r") as f:
            st.session_state.settings = json.load(f)

if __name__ == '__main__':
    local = not st.session_state.settings["online_deployment"]
    # Sections follow the path a user takes: set up and run, read the main
    # answers, dig deeper. Expert and quality-control pages come last.
    pages = {
        "Start here": [
            st.Page(Path("content", "quickstart.py"), title="Home", icon="🏠", default=True),
        ],
        "Step 1 · Set up and run": [
            st.Page(Path("content", "workflow_fileupload.py"), title="Upload data", icon="📁"),
            st.Page(Path("content", "workflow_configure.py"), title="Describe experiment", icon="⚙️"),
            st.Page(Path("content", "workflow_run.py"), title="Run analysis", icon="🚀"),
        ],
        "Step 2 · Main results": [
            st.Page(Path("content", "results_overview.py"), title="Run summary", icon="🧾"),
            st.Page(Path("content", "compare_groups.py"), title="Compare groups", icon="⚖️"),
            st.Page(Path("content", "results_volcano.py"), title="Volcano plot", icon="🌋"),
            st.Page(Path("content", "results_proteins.py"), title="Identified proteins", icon="🧩"),
        ],
        # The sidebar shows the first ten pages and folds the rest behind
        # "View more", so everything from Clustered heatmap on is one click away.
        "Step 3 · Explore further": [
            st.Page(Path("content", "results_pca.py"), title="Sample similarity (PCA)", icon="📊"),
            st.Page(Path("content", "enrichment.py"), title="GO enrichment", icon="📉"),
            st.Page(Path("content", "results_heatmap_clustered.py"), title="Clustered heatmap", icon="🧬"),
            st.Page(Path("content", "results_abundance.py"), title="Protein abundances", icon="📋"),
            st.Page(Path("content", "results_samples.py"), title="Per-sample counts", icon="🧪"),
            st.Page(Path("content", "results_heatmap.py"), title="Heatmap", icon="🔥"),
        ],
        "Export": [
            # A local install exports results to a folder instead of downloading them.
            st.Page(Path("content", "results_download.py"), title="Runs and exports" if local else "Runs and downloads", icon="🗂️"),
        ],
        # The four steps Compare groups runs with recommended settings, one at a time.
        "Fine-tune the comparison": [
            st.Page(Path("content", "filtering.py"), title="1. Filter proteins", icon="🧹"),
            st.Page(Path("content", "imputation.py"), title="2. Missing values", icon="🩹"),
            st.Page(Path("content", "normalization.py"), title="3. Normalize", icon="⚖️"),
            st.Page(Path("content", "statistical.py"), title="4. Statistical test", icon="🔢"),
        ],
        "Quality control and details": [
            st.Page(Path("content", "results_qc.py"), title="Identification quality", icon="🩺"),
            st.Page(Path("content", "results_filtered.py"), title="Peptide-spectrum matches", icon="🎯"),
            st.Page(Path("content", "results_mass_shifts.py"), title="Mass shifts", icon="📐"),
            st.Page(Path("content", "results_rescoring.py"), title="Rescored hits (Percolator)", icon="📈"),
            st.Page(Path("content", "results_database_search.py"), title="Raw search hits (Comet)", icon="🔬"),
        ],
    }

    pg = st.navigation(pages)
    pg.run()

