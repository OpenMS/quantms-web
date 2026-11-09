# Quality control coverage

Which quality control metrics of a DDA run the Viewer shows, and which it does not show yet.
Kept up to date as metrics are added; the open ones are listed roughly by value.

Status: **shown**, **partly** (shown, but not per file or not in every mode), **missing**.

## Identification

| Metric | Status | Where |
|---|---|---|
| Target and decoy score distribution, FDR cutoff | shown | Identification quality |
| PSMs per file | shown | Identification quality, Run summary |
| ID rate per file (PSMs per MS2 spectrum) and files that fall short | shown | Run summary, Diagnostics |
| PSMs, peptides and proteins per sample | shown | Per-sample counts |
| Precursor mass error | partly | Identification quality (all files together; not per file, not over retention time) |
| Precursor mass shifts and the modifications they match | shown | Mass shifts |
| Precursor charge states per file | shown | Charge states and missed cleavages |
| Missed cleavages per file | shown | Charge states and missed cleavages |
| Contaminants: share of PSMs per file | partly | Diagnostics (a cause); no share of intensity per file |
| Fragment mass error | missing | |
| Peptide length distribution | missing | |
| Identifications over retention time (per file) | missing | |
| MS2 spectra identified more than once (oversampling) | missing | |
| Modified peptides per file (oxidation, others searched) | missing | |

## Acquisition

| Metric | Status | Where |
|---|---|---|
| MS1 and MS2 spectra per file | partly | Run summary (MS2 counts) |
| Total ion current over retention time, per file | missing | |
| MS2 spectra over retention time (duty cycle) | missing | |
| Ion injection time of MS1 and MS2 spectra | missing | |
| Precursor intensity of identified and unidentified MS2 spectra | missing | |
| Chromatographic peak width of identified peptides | missing | |

## Quantification

| Metric | Status | Where |
|---|---|---|
| Intensity distribution per sample, before and after normalization | shown | Normalization |
| Missing values per sample and by intensity | shown | Missing values |
| Sample similarity (PCA, clustered heatmap) | shown | PCA, Clustered heatmap |
| Retention time alignment shift per file | missing | |
| Features transferred between runs (match between runs) per file | missing | |
| Contaminant share of the total intensity per sample | missing | |
| Reporter ion intensity per channel (TMT/iTRAQ) | partly | Normalization shows the channel distributions; no per-file reporter view |
| Missing reporter ions per channel | missing | |
| Isolation purity of labelled precursors | missing | |

## Next candidates

1. Identifications and precursor mass error over retention time, per file (calibration drift, gradient problems).
2. Total ion current and MS2 spectra over retention time, per file.
3. Contaminant share of the intensity per sample.
4. Reporter ion intensity and missing reporter ions per channel for TMT/iTRAQ runs.
