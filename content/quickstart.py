"""
OpenDDA Quickstart Page.

This page provides an overview of the DDA quantification workflows (label-free
and TMT) and guidance for getting started with the analysis pipeline.
"""

from pathlib import Path

import streamlit as st
from src.common.common import page_setup

WINDOWS_APP_PATH = Path("/app/OpenMS-App.zip")


@st.cache_resource
def load_windows_app_bytes() -> bytes | None:
    if WINDOWS_APP_PATH.exists():
        return WINDOWS_APP_PATH.read_bytes()
    return None


def render_windows_download_box(app_bytes: bytes) -> None:
    """Render a styled offline-download card for the Windows app."""
    container_key = "windows_download_container"

    st.markdown(
        """
        <h4 style="color: #6c757d; margin-bottom: 1rem; font-size: 1.3rem; font-weight: 600; text-align: center;">
            Want to run free and open DDA analysis (LFQ and TMT) offline?
        </h4>
        """,
        unsafe_allow_html=True,
    )

    with st.container(key=container_key):
        st.markdown(
            """
            <h4 style="color: #6c757d; margin-bottom: 0.75rem; font-size: 1.1rem; font-weight: 600;">
                OpenDDA for Windows
            </h4>
            <p style="color: #6c757d; margin-bottom: 1rem;">
                You can download an offline version for Windows systems below.
            </p>
            """,
            unsafe_allow_html=True,
        )

        cols = st.columns([2, 3, 2])
        with cols[1]:
            st.download_button(
                label="📥 Download for Windows",
                data=app_bytes,
                file_name="OpenMS-App.zip",
                mime="application/zip",
                type="secondary",
                use_container_width=True,
                help="Download OpenDDA for Windows systems",
            )

        st.markdown(
            """
            <div style="text-align: center; margin-top: 1rem; color: #6c757d;">
                Extract the zip file and run the installer (.msi) to install the app.
                Launch using the desktop icon after installation.<br>
                Even offline, it's still a web app - just packaged so you can use it without an internet connection.
            </div>
            """,
            unsafe_allow_html=True,
        )

    st.markdown(
        f"""
        <style>
        .st-key-{container_key} {{
            background: linear-gradient(135deg, #f8f9fa 0%, #f1f3f4 100%) !important;
            border: 1px solid #e0e0e0 !important;
            border-radius: 8px !important;
            padding: 1.5rem !important;
            margin: 1rem 0 !important;
            text-align: center !important;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05) !important;
        }}

        .st-key-{container_key} > div {{
            background: transparent !important;
        }}
        </style>
        """,
        unsafe_allow_html=True,
    )


page_setup(page="main")

st.markdown("# OpenDDA: DDA Quantitative Proteomics")

windows_app_bytes = load_windows_app_bytes()
if windows_app_bytes is not None:
    render_windows_download_box(windows_app_bytes)

st.markdown(
    """
This application provides complete **Data-Dependent Acquisition (DDA)** quantification
workflows for proteomics data analysis. The pipeline identifies and quantifies proteins
from mass spectrometry data and supports two quantification strategies:

- **Label-free quantification (LFQ)**: compares precursor intensities across separately measured samples.
- **Tandem Mass Tag (TMT) quantification**: compares reporter-ion intensities of multiplexed samples measured in one run.

Choose the **Analysis Mode** (LFQ or TMT) on the Configure page.
"""
)

st.info(
    "These workflows mirror the **DDA-LFQ and DDA-ISO (TMT) branches of the quantms Nextflow workflow**."
)

st.markdown("## Workflow Overview")

lfq_col, tmt_col = st.columns(2)

with lfq_col:
    st.markdown(
        """
### Label-free (LFQ)

| Stage | Tool |
|-------|------|
| **1. Identification** | Comet |
| **2. Rescoring** | Percolator |
| **3. Filtering** | IDFilter |
| **4. Quantification** | ProteomicsLFQ |
| **5. Statistical Analysis** | Built-in |
"""
    )

with tmt_col:
    st.markdown(
        """
### TMT

| Stage | Tool |
|-------|------|
| **1. Reporter extraction** | IsobaricAnalyzer |
| **2. Identification** | Comet |
| **3. Rescoring** | Percolator |
| **4. Filtering** | IDFilter |
| **5. Protein inference** | ProteinInference |
| **6. Quantification** | ProteinQuantifier |
| **7. Statistical Analysis** | Built-in |
"""
    )

st.markdown("## Getting Started")

st.markdown(
    """
Follow these steps to run your analysis:

### 1. Upload Files
Upload your mzML mass spectrometry files and a protein FASTA database.
"""
)
st.page_link("content/workflow_fileupload.py", label="Go to File Upload", icon="📁")

st.markdown(
    """
### 2. Configure Parameters
Choose LFQ or TMT, then set up search parameters, sample groups (or TMT channels), and analysis settings.
"""
)
st.page_link("content/workflow_configure.py", label="Go to Configure", icon="⚙️")

st.markdown(
    """
### 3. Run Workflow
Execute the analysis pipeline and monitor progress.
"""
)
st.page_link("content/workflow_run.py", label="Go to Run", icon="🚀")

st.markdown(
    """
### 4. Explore Results
View identification results, quantification tables, and statistical visualizations.
- **Database Search**: View Comet PSM identification results
- **Rescoring**: Examine Percolator statistical validation output
- **Filtered PSMs**: Inspect FDR-controlled peptide identifications
- **Abundance**: Protein and PSM quantification tables with statistics
- **Volcano Plot**: Differential expression analysis visualization
- **PCA**: Principal component analysis of sample relationships
- **Heatmap**: Hierarchically clustered expression patterns
- **Spectral Library**: Download the generated spectral library for DIA/SWATH analysis
"""
)

st.markdown("## Workspaces")

st.markdown(
    """
**Workspaces** store your inputs, parameters, and results for each analysis session.

- In **online mode**: The workspace ID is embedded in the URL for easy sharing
- In **local mode**: Create and manage separate workspaces for different projects

Your workspace persists between sessions, allowing you to return to a running or completed
workflow at any time.
"""
)
