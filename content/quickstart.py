"""
OpenDDA Quickstart Page.

This page provides an overview of the DDA quantification workflows (label-free
and TMT) and guidance for getting started with the analysis pipeline.
"""

from pathlib import Path

import streamlit as st
from src.common.common import page_setup
from src.common.guide import page_link

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

st.markdown("# OpenDDA: find out which proteins change in your samples")

windows_app_bytes = load_windows_app_bytes()
if windows_app_bytes is not None:
    render_windows_download_box(windows_app_bytes)

st.markdown(
    """
OpenDDA turns raw mass spectrometry files from a **data-dependent acquisition (DDA)**
experiment into a list of proteins, how much of each protein every sample contains,
and which proteins differ between your experimental groups. You do not need to know
the search engines behind it: every setting starts at a value that suits a typical
experiment, and each page explains what you are looking at.
"""
)

st.markdown("## How it works")
steps = st.columns(3)
with steps[0]:
    st.markdown(
        """
### 1 · Set up and run
Upload your **mzML** files and a **FASTA** file with the protein sequences of your
organism (or import a public dataset from PRIDE). Say which sample belongs to which
group, then press Run. A run takes minutes to hours depending on the data.
"""
    )
    page_link("upload", "Upload data")
with steps[1]:
    st.markdown(
        """
### 2 · Read the main results
The **run summary** tells you whether the run went well. **Compare groups** and the
**volcano plot** show the proteins that go up or down between your groups.
**Identified proteins** lists everything that was found.
"""
    )
    page_link("compare", "Compare groups")
with steps[2]:
    st.markdown(
        """
### 3 · Explore further
Check that replicates look alike (**PCA**, **clustered heatmap**), and find biological
processes over-represented among the changed proteins (**GO enrichment**).
"""
    )
    page_link("enrichment", "GO enrichment")

st.markdown("## Which page answers my question?")
st.markdown(
    """
| Your question | Page |
|---|---|
| Did the run work, and how many proteins were found? | **Run summary** |
| Which proteins are up or down between my groups? | **Compare groups**, **Volcano plot** |
| Is my protein of interest there, and which peptides support it? | **Identified proteins** |
| How much of each protein is in each sample? | **Protein abundances** |
| Do my replicates agree? Is there an outlier sample? | **Sample similarity (PCA)**, **Clustered heatmap** |
| What biology do the changed proteins point to? | **GO enrichment** |
| I want the tables for my own analysis. | **Download results** |
"""
)
st.caption(
    "The sections **Fine-tune the comparison** and **Quality control and details** in the sidebar "
    "are for adjusting the statistics step by step and for checking the identifications in depth. "
    "You can ignore them for a first analysis."
)

st.markdown("## Two kinds of experiments")
lfq_col, tmt_col = st.columns(2)
with lfq_col:
    st.markdown(
        """
### Label-free (LFQ)
Each sample was measured in its own run. Proteins are compared by the intensity of
their peptides across runs. Choose this if you did not label your samples.
"""
    )
with tmt_col:
    st.markdown(
        """
### Tandem Mass Tags (TMT)
Several samples were labelled with isobaric tags, pooled and measured together.
Proteins are compared by the tag (reporter ion) intensities within each run.
"""
    )
st.caption("You choose the experiment type on the Describe experiment page.")

with st.expander("📖 What the analysis does behind the scenes"):
    st.markdown(
        """
The pipeline mirrors the DDA-LFQ and DDA-ISO (TMT) branches of the
[quantms](https://quantms.org) Nextflow workflow and uses OpenMS tools throughout:

| Stage | What it does | Tool |
|---|---|---|
| Identification | Matches each fragmentation spectrum to a peptide from your FASTA file | Comet |
| Rescoring | Re-ranks the matches with machine learning to separate true from false ones | Percolator |
| Filtering | Keeps matches at a 1 % false discovery rate, judged against reversed "decoy" sequences | IDFilter |
| Quantification (LFQ) | Measures each peptide's signal in every run and sums it per protein | ProteomicsLFQ |
| Quantification (TMT) | Reads the reporter-ion intensities and sums them per protein | IsobaricAnalyzer, ProteinQuantifier |
| Group comparison | Filters, fills in missing values, normalizes and tests every protein | Built in |
"""
    )

with st.expander("💾 Workspaces: where your data and results live"):
    st.markdown(
        """
A **workspace** keeps your files, settings and results together.

- **Online**, the workspace ID is part of the page address. Bookmark it to come back
  to a running or finished analysis, or share it with a colleague.
- **On your own computer**, create a separate workspace per project in the sidebar.
"""
    )
