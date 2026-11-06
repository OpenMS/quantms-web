"""
OpenDDA Quickstart Page.

This page provides an overview of the DDA quantification workflows (label-free
and TMT) and guidance for getting started with the analysis pipeline.
"""

import streamlit as st
from src.common.common import get_branding, page_setup
from src.common.guide import page_link
from src.common.installers import installers_available, render_installer_buttons

def render_desktop_download_box() -> None:
    """Render the offline-download card: one button per platform with an installer."""
    if not installers_available():
        return
    container_key = "desktop_download_container"

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
                OpenDDA desktop app
            </h4>
            <p style="color: #6c757d; margin-bottom: 1rem;">
                Download the offline version for Windows, macOS or Linux below.
            </p>
            """,
            unsafe_allow_html=True,
        )

        render_installer_buttons(app_name="OpenDDA")

    st.markdown(
        f"""
        <style>
        .st-key-{container_key} {{
            background: rgba(128, 128, 128, 0.08) !important;
            border: 1px solid rgba(128, 128, 128, 0.3) !important;
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

title_col, logo_col = st.columns([6, 1], vertical_alignment="center")
with title_col:
    st.markdown("# OpenDDA: Find differentially abundant proteins in your samples")
with logo_col:
    st.image(get_branding()["app_icon"], width=120)

render_desktop_download_box()

st.markdown(
    """
OpenDDA turns raw mass spectrometry files from a **data-dependent acquisition (DDA)**
experiment into a list of proteins, how much of each protein every sample contains,
and which proteins are differentially abundant between your experimental groups. You do not need to know
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
Choose your runs (**mzML**, Thermo *.raw* or Bruker *.d*) and the protein sequences of
your organism (a **FASTA** file or a UniProt proteome), or reanalyse a public dataset
from PRIDE. Say which sample belongs to which group, then press Run. A run takes minutes
to hours depending on the data.
"""
    )
    page_link("configure", "Describe experiment")
with steps[1]:
    st.markdown(
        """
### 2 · Read the main results
The **Viewer** opens on the **run summary**, which tells you whether the run went well,
the **group comparison** and **volcano plot**, which show the proteins that go up or
down between your groups, and the **identified proteins**.
"""
    )
    page_link("viewer", "Viewer")
with steps[2]:
    st.markdown(
        """
### 3 · Explore further
**Edit layout** in the Viewer adds more components: check that replicates look alike
(**PCA**, **clustered heatmap**), find biological processes over-represented among the
changed proteins (**GO enrichment**), or put two runs side by side.
"""
    )
    page_link("viewer", "Viewer")

st.markdown("## Which component answers my question?")
st.markdown(
    """
| Your question | Viewer component |
|---|---|
| Did the run work, and how many proteins were found? | **Run summary** |
| Which proteins are up or down between my groups? | **Group comparison**, **Volcano plot** |
| Is my protein of interest there, and which peptides support it? | **Identified proteins** |
| How much of each protein is in each sample? | **Protein and peptide abundances** |
| Do my replicates agree? Is there an outlier sample? | **Sample similarity (PCA)**, **Clustered heatmap** |
| What biology do the changed proteins point to? | **GO enrichment** |
| I want the tables for my own analysis. | **Export**, below the Viewer |
"""
)
st.caption(
    "**Edit comparison** adjusts the statistics step by step, and the components under "
    "**Identification** (quality, peptide-spectrum matches, mass shifts, rescored and raw hits) check "
    "the identifications in depth. You can ignore them for a first analysis."
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
| Identification | Matches each fragmentation spectrum to a peptide from your FASTA file | ProSE |
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
