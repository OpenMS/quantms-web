# OpenDDA (quantms-web)

OpenDDA's own guidance, imported by the template-owned `CLAUDE.md` (which describes the framework).

**OpenDDA** is quantitative DDA proteomics: label-free (LFQ) and isobaric (TMT/iTRAQ) quantification. The workflow (`src/WorkflowTest.py`, class `WorkflowTest`) runs database search, Percolator rescoring, FDR filtering and ProteomicsLFQ (or the isobaric path), followed by statistics: normalization, imputation, group comparison (volcano), PCA, heatmaps, GO and pathway enrichment.

- Pages: `content/workflow_*.py` (upload, configure, run), `content/results_*.py` (one per results view), and the statistics pages (`normalization.py`, `imputation.py`, `statistical.py`, `compare_groups.py`, `enrichment.py`).
- Results logic lives in `src/common/` (`protein_table.py`, `spectrum_view.py`, `qc_plots.py`, `samples.py`, `differential.py`, `isobaric.py`, `mass_shifts.py`, `contaminants.py`, `ms2_prediction.py`, ...). Tables go through `show_insight_table(df, name=...)` from `src/common/insight_tables.py`, never `st.dataframe()` on full-size frames.
- PRIDE import (accessions, RAW conversion, SDRF) is `src/pride.py` / `src/pride_ui.py`, next to the template's Import from a link (`src/data_import/`).
- Insight components request a rerun on first render, so never draw them only inside `if st.button(...)`: store the result in `st.session_state` and render outside the button block.
- Deployed as `quantms-ddalfq` (`k8s/overlays/prod`), hosts `opendda.webapps.openms.de` / `.org`. Release tags look like `OpenDDA-1.3.0`.
