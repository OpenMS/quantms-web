import streamlit as st
from src.common.common import page_setup
from src.common.guide import page_intro
from src.WorkflowTest import WorkflowTest


params = page_setup()

wf = WorkflowTest()

st.title("Describe your experiment")
page_intro(
    """
For most experiments only three things need your attention:

1. **Instrument preset** (buttons below): sets the mass tolerances for your mass spectrometer.
2. **Experiment type**: label-free (LFQ) or labelled (TMT/iTRAQ).
3. **Samples and groups** (first tab): which sample belongs to which condition.

Everything else starts at values that suit a typical experiment and is saved as you go.
Switch on **Show advanced parameters** only if you need to change a specific setting.
""",
    background="""
- **mzML** is an open file format for raw mass spectrometry data. Convert vendor files
  (e.g. Thermo *.raw*) with [ProteoWizard msConvert](https://proteowizard.sourceforge.io/).
- The **FASTA** file lists the protein sequences of your organism, e.g. the reference
  proteome from [UniProt](https://www.uniprot.org/proteomes). Only proteins in this file can be found.
- **Modifications** like phosphorylation are added under *Database search*. Each extra
  variable modification makes the search slower, so add only those you are interested in.
""",
    background_title="Input files and terms explained",
)

wf.show_parameter_section()
