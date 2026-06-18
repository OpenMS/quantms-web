# quantms-web — DDA Label-Free Quantification

A browser-based **Data-Dependent Acquisition (DDA) Label-Free Quantification** workflow for proteomics. Upload mzML files and a protein FASTA; get identified and quantified proteins with volcano plots, PCA, and clustered heatmaps. No CLI, no Nextflow config.

quantms-web mirrors the **dda-lfq branch of the [quantms Nextflow workflow](https://github.com/bigbio/quantms)** but runs as a [Streamlit](https://streamlit.io) app powered by OpenMS TOPP tools.

## Pipeline

| Stage | Tool | What it does |
|---|---|---|
| 1. Identification | Comet | Peptide-spectrum matching against a protein database |
| 2. Rescoring | Percolator | ML-based statistical validation of PSMs |
| 3. Filtering | IDFilter | FDR-controlled peptide identification filtering |
| 4. Quantification | ProteomicsLFQ | Label-free quantification across samples |
| 5. Analysis | Built-in | Volcano plots, PCA, heatmaps, spectral library export |

## Run locally

Install the Python dependencies and launch:

```bash
git clone https://github.com/OpenMS/quantms-web.git
cd quantms-web
pip install -r requirements.txt
streamlit run app.py
```

The full pipeline runs locally once the [OpenMS Command Line Tools](https://openms.readthedocs.io/en/latest/about/installation.html) are on your `PATH` — they provide Comet, Percolator, ProteomicsLFQ, and the rest of the TOPP suite. With Python alone, the pyOpenMS-backed parts of the UI still work.

## Run with Docker

Ships OpenMS and the search engines together so the full pipeline works out of the box:

```bash
docker-compose up -d --build
```

Open http://localhost:8501.

## Windows installer

Download the latest `.msi` from [Releases](https://github.com/OpenMS/quantms-web/releases) and double-click to install. Standalone — no Python or Docker required.

## Workspaces

Every analysis session runs in an isolated **workspace** that persists inputs, parameters, and results. In online deployments the workspace ID is part of the URL, so runs are resumable and shareable.

## Citation

Müller, T. D., Siraj, A., et al. *OpenMS WebApps: Building User-Friendly Solutions for MS Analysis.* Journal of Proteome Research (2025). [doi:10.1021/acs.jproteome.4c00872](https://doi.org/10.1021/acs.jproteome.4c00872)

## References

- Pfeuffer, J., Bielow, C., Wein, S. et al. *OpenMS 3 enables reproducible analysis of large-scale mass spectrometry data.* Nat Methods 21, 365–367 (2024). [doi:10.1038/s41592-024-02197-7](https://doi.org/10.1038/s41592-024-02197-7)
- Röst HL, Schmitt U, Aebersold R, Malmström L. *pyOpenMS: a Python-based interface to the OpenMS mass-spectrometry algorithm library.* Proteomics 14, 74–77 (2014). [doi:10.1002/pmic.201300246](https://doi.org/10.1002/pmic.201300246)
