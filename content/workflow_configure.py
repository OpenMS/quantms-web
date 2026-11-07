import streamlit as st
from src.common.common import page_setup
from src.common.guide import PAGES, page_intro
from src.WorkflowTest import WorkflowTest


params = page_setup()

wf = WorkflowTest()

st.title("Describe your experiment")
page_intro(
    """
Choose your runs and your analysis in the box below, then say which sample belongs to which
group in **Samples and groups**. All other settings start at values that suit a typical experiment.

- **MS data**: one file per run, mzML, Thermo *.raw* or a Bruker timsTOF *.d* folder
  (vendor files are converted for you). **📂 Browse** opens your files by folder and adds new ones.
  For a TMT experiment, choose every run (fraction) of the labelled pool.
- **Protein database**: a FASTA file, or a UniProt proteome chosen with **📂 Browse**. Decoys
  and common contaminants are added for you.
- **Instrument**: sets the mass tolerances of the search for your mass spectrometer. The
  instrument the runs were measured on is read from the files.

No data at hand? A small demo dataset (yeast with spiked-in human proteins in two
concentrations) is used until you add your own.
""",
    background="""
- **mzML** is an open file format for raw mass spectrometry data. Formats other than Thermo
  *.raw* and Bruker *.d* convert to it with [ProteoWizard msConvert](https://proteowizard.sourceforge.io/)
  (keep *peak picking* on).
- The **FASTA** file lists the protein sequences of your organism, e.g. the reference
  proteome from [UniProt](https://www.uniprot.org/proteomes). Only proteins in this file can be found.
- **Modifications** like phosphorylation are added under *Database search*. Each extra
  variable modification makes the search slower, so add only those you are interested in.
- An **open search** accepts any precursor mass within ±500 Da of a peptide, so peptides with
  modifications nobody searched for are found too.
""",
    background_title="Input files and terms explained",
)

wf.analysis_section()
wf.show_parameter_section()

# A bar at the bottom of the window: opens Run analysis and starts the analysis there
wf.show_run_button(PAGES["run"][0], label="Run analysis", summary=wf.run_bar_summary(),
                   help="Save these settings, open Run analysis and start the analysis of the selected runs.")
