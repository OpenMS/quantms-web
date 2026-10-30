import streamlit as st
from src.common.common import page_setup
from src.common.guide import page_intro
from src.WorkflowTest import WorkflowTest

params = page_setup()
wf = WorkflowTest()

st.title("Upload your data")
page_intro(
    """
You need two kinds of files:

- **MS data (mzML)**: one file per measured run. Convert vendor raw files to mzML first.
  For a TMT experiment, upload every run (fraction) of the labelled pool; the groups are
  then set per label channel, not per file.
- **Protein database (FASTA)**: the protein sequences of your organism. Decoys and common
  contaminants are added automatically.

No data at hand? The app comes with a small demo dataset (yeast with spiked-in human
proteins in two concentrations) that is used when nothing is uploaded. You can also
**import a public dataset from PRIDE** to reanalyse it.
""",
    background="""
- Convert Thermo *.raw*, Bruker *.d* and other vendor formats to mzML with
  [ProteoWizard msConvert](https://proteowizard.sourceforge.io/) (keep *peak picking* on).
- Download a reference proteome as FASTA from [UniProt](https://www.uniprot.org/proteomes),
  ideally *one protein sequence per gene* to keep the search fast.
""",
    background_title="How to get these files",
)

wf.show_file_upload_section()
