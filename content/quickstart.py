"""Quickstart page for the DDA LFQ and TMT proteomics workflows."""

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
            Want to run free and open DDA proteomics analysis offline?
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

st.markdown("# DDA Quantitative Proteomics: LFQ & TMT")

windows_app_bytes = load_windows_app_bytes()
if windows_app_bytes is not None:
    render_windows_download_box(windows_app_bytes)

st.markdown(
    """
This application supports two complete **Data-Dependent Acquisition (DDA) quantitative
proteomics workflows**:

- **LFQ (Label-Free Quantification)** for comparing protein abundance across independently
  acquired samples without isobaric labels.
- **TMT (Tandem Mass Tag)** for reporter-ion-based multiplexed quantification across channels
  within an isobarically labeled experiment.

Both workflows start from **mzML files and a protein FASTA database**, perform peptide and
protein identification and quantification, and feed a shared downstream analysis workflow for
filtering, missing-value imputation, normalization, statistical inference, and visualization.
"""
)

st.info(
    "Select **LFQ** or **TMT** from Analysis Mode on the Configure page. "
    "The available parameters and execution pipeline change automatically for the selected mode."
)

st.markdown("## LFQ Workflow")

st.markdown(
    """
Use LFQ when each biological sample was acquired as a separate, unlabeled LC-MS/MS run.

| Stage | Tool | Description |
|-------|------|-------------|
| **1. Decoy preparation** | DecoyDatabase | Generates a target-decoy FASTA when the uploaded database does not already contain decoys |
| **2. Identification** | CometAdapter | Searches each mzML file against the protein database |
| **3. Rescoring** | PercolatorAdapter | Improves PSM confidence and assigns statistical scores |
| **4. Identification filtering** | IDFilter | Retains identifications according to the configured score/FDR threshold |
| **5. Spectral library (optional)** | EasyPQP | Builds a spectral library for downstream DIA/SWATH use |
| **6. Label-free quantification** | ProteomicsLFQ | Aligns and quantifies peptide/protein signals across mzML runs |

On the **Group Selection** tab, assign each mzML file to a biological group such as
`control` or `case`. These assignments are used by the downstream statistical pages.
"""
)

st.markdown("## TMT Workflow")

st.markdown(
    """
Use TMT when samples were labeled with tandem mass tags and combined into a multiplexed run.

| Stage | Tool | Description |
|-------|------|-------------|
| **1. Reporter-ion extraction** | IsobaricAnalyzer | Extracts and quantifies reporter-ion signals for the selected TMT plex |
| **2. Identification** | CometAdapter | Searches MS/MS spectra against the protein database with TMT modifications |
| **3. Rescoring** | PercolatorAdapter | Rescores peptide-spectrum matches and estimates confidence |
| **4. PSM filtering** | IDFilter | Applies the configured identification threshold |
| **5. ID mapping** | IDMapper | Maps filtered peptide identifications to isobaric consensus features |
| **6. Run merging** | FileMerger | Combines mapped consensus results across input runs |
| **7. Protein inference** | ProteinInference | Infers protein groups from identified peptides |
| **8. Protein filtering** | IDFilter | Applies protein/peptide-level filtering and removes decoys |
| **9. Conflict resolution** | IDConflictResolver | Resolves conflicting peptide-to-protein assignments |
| **10. Protein quantification** | ProteinQuantifier | Produces the final channel-level protein abundance matrix |

Choose the isobaric type in **IsobaricAnalyzer**, then assign each TMT channel to a biological
group on the **Group Selection** tab. Enter `skip` for an unused channel.
"""
)

st.markdown("## Getting Started")

st.markdown(
    """
Follow these steps to run your analysis:

### 1. Upload Files
Upload one or more mzML mass spectrometry files and a protein FASTA database. The same upload
page is used for LFQ and TMT analyses.
"""
)
st.page_link("content/workflow_fileupload.py", label="Go to File Upload", icon="📁")

st.markdown(
    """
### 2. Configure Parameters
Select **LFQ** or **TMT** as the Analysis Mode. Configure the tools shown for that workflow,
then assign mzML files (LFQ) or reporter channels (TMT) to biological groups.
"""
)
st.page_link("content/workflow_configure.py", label="Go to Configure", icon="⚙️")

st.markdown(
    """
### 3. Run Workflow
Execute the selected workflow and monitor each processing stage. Results are written to the
current workspace so they remain available after the run finishes.
"""
)
st.page_link("content/workflow_run.py", label="Go to Run", icon="🚀")

st.markdown("### 4. Review Identification and Quantification Results")

st.markdown(
    """
The **Results** section exposes the main outputs from either workflow:

- **Database Search**: inspect Comet peptide-spectrum matches, peptide sequences, and MS2 spectra.
- **Rescoring**: review Percolator-rescored identifications.
- **Filtered PSMs**: inspect identifications retained after score/FDR filtering.
- **Abundance**: review the final protein abundance matrix. LFQ displays run-level intensities;
  TMT displays reporter-channel-level protein intensities.
"""
)

st.page_link("content/results_abundance.py", label="Go to Abundance", icon="📊")

st.markdown("### 5. Run Differential Protein Analysis")

st.markdown(
    """
After quantification, continue through the downstream analysis pages in order:

| Step | Available analysis |
|------|--------------------|
| **Filtering** | Remove proteins with low abundance, low repeatability, or low variance |
| **Imputation** | Handle missing values using MAR mean/median strategies or MNAR small-value replacement |
| **Normalization** | Apply mathematical transformation, sample normalization, and optional row scaling |
| **Statistical** | Test group differences and apply BH, Bonferroni, or no multiple-testing correction |
| **Volcano** | Explore effect size and adjusted-significance thresholds |
| **PCA** | Inspect sample separation using high-variance proteins |
| **Heatmap** | View abundance patterns across proteins and samples |
| **Clustered Heatmap** | Explore hierarchically clustered proteins and samples |
| **Pathway Analysis** | Run GO enrichment on statistically significant proteins |

Each preprocessing page uses the output of the preceding step when available. The recommended
order is **Filtering → Imputation → Normalization → Statistical → Visualization/Pathway Analysis**.
"""
)

analysis_links = st.columns(4)
with analysis_links[0]:
    st.page_link("content/filtering.py", label="Filtering", icon="🔍")
with analysis_links[1]:
    st.page_link("content/imputation.py", label="Imputation", icon="🧩")
with analysis_links[2]:
    st.page_link("content/normalization.py", label="Normalization", icon="📐")
with analysis_links[3]:
    st.page_link("content/statistical.py", label="Statistical", icon="🧪")

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
