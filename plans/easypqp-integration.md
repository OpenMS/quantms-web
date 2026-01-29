# EasyPQP Spectral Library Integration

## Summary

Integrate EasyPQP into the workflow to generate OpenSWATH-compatible spectral libraries (TSV format) from filtered PSMs. Adds configuration options, execution steps after IDFilter, and download button on results page.

---

## Requirements (Confirmed with User)

| Requirement | Decision |
|-------------|----------|
| Input source | IDFilter results (filtered idXML) + source mzML |
| Output format | TSV for OpenSWATH |
| Merging | Merge all samples into single library |
| FDR | Configurable parameter (new tab) |
| Decoys | Checkbox + method dropdown |
| UI location | Config tab + download on results_filtered.py |
| Execution timing | After IDFilter, before ProteomicsLFQ |

---

## Files to Modify

| File | Changes |
|------|---------|
| `src/WorkflowTest.py` | Add config tab (line ~45), add execution step (after line 634) |
| `content/results_filtered.py` | Add download section (before line 75) |
| `default-parameters.json` | Add default values for new parameters |

---

## 1. Configuration Changes (`src/WorkflowTest.py`)

### 1.1 Add New Tab (line 45)

**Before:**
```python
t = st.tabs(["**Identification**", "**Rescoring**", "**Filtering**", "**Quantification**", "**Group Selection**"])
```

**After:**
```python
t = st.tabs(["**Identification**", "**Rescoring**", "**Filtering**", "**Library Generation**", "**Quantification**", "**Group Selection**"])
```

**Note:** Tab indices shift: t[3]=Library, t[4]=Quantification, t[5]=Group Selection

### 1.2 Add Tab Content (after line 165, before current `with t[3]`)

Insert new `with t[3]:` block for Library Generation:

```python
with t[3]:  # Library Generation
    st.info("""
    **Spectral Library Generation (EasyPQP):**
    Generate a spectral library from filtered PSMs for targeted proteomics (DIA/SWATH).
    The library is created from all filtered idXML files combined with their source mzML spectra.
    """)

    self.ui.input_widget(
        key="generate-library",
        default=False,
        name="Generate Spectral Library",
        widget_type="checkbox",
        help="Enable spectral library generation using EasyPQP.",
        reactive=True,
    )

    self.params = self.parameter_manager.get_parameters_from_json()

    if self.params.get("generate-library", False):
        st.markdown("---")
        st.markdown("**FDR Options**")

        self.ui.input_widget(
            key="library-use-fdr",
            default=False,
            name="Apply additional FDR filtering in EasyPQP",
            widget_type="checkbox",
            help="If disabled (recommended), uses --nofdr since IDFilter already applied FDR control.",
            reactive=True,
        )

        self.params = self.parameter_manager.get_parameters_from_json()

        if self.params.get("library-use-fdr", False):
            self.ui.input_widget(
                key="library-psm-fdr",
                default=0.01,
                name="PSM FDR Threshold",
                widget_type="number",
                min_value=0.001,
                max_value=1.0,
                step_size=0.01,
                help="FDR threshold for PSMs (e.g., 0.01 = 1% FDR).",
            )

        st.markdown("---")
        st.markdown("**Decoy Generation**")

        self.ui.input_widget(
            key="library-generate-decoys",
            default=True,
            name="Generate Decoy Library",
            widget_type="checkbox",
            help="Generate decoy transitions for FDR control in OpenSWATH.",
            reactive=True,
        )

        self.params = self.parameter_manager.get_parameters_from_json()

        if self.params.get("library-generate-decoys", True):
            self.ui.input_widget(
                key="library-decoy-method",
                default="shuffle",
                name="Decoy Method",
                widget_type="selectbox",
                options=["shuffle", "reverse", "pseudo-reverse", "shift"],
                help="Method for generating decoy sequences.",
            )
```

### 1.3 Update Tab Index References

Update all subsequent tab references:
- `with t[3]:` (Quantification) becomes `with t[4]:`
- `with t[4]:` (Group Selection) becomes `with t[5]:`

---

## 2. Execution Changes (`src/WorkflowTest.py`)

### Insert After Line 634 (after `self.logger.log("Filtering complete")`)

```python
# ================================
# 3.5 EasyPQP Spectral Library Generation (optional)
# ================================
if self.params.get("generate-library", False):
    self.logger.log("Building spectral library with EasyPQP...")
    st.info("Generating spectral library from filtered PSMs...")

    library_dir = Path(self.workflow_dir, "results", "library")
    library_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Convert each filtered idXML + mzML to psms/peaks
    psms_files = []
    peaks_files = []

    for filter_idxml in filter_results:
        filter_stem = Path(filter_idxml).stem  # e.g., "sample1_filter"
        original_stem = filter_stem.replace("_filter", "")

        # Find matching mzML
        matching_mzml = None
        for mzml in in_mzML:
            if Path(mzml).stem == original_stem:
                matching_mzml = mzml
                break

        if matching_mzml is None:
            self.logger.log(f"WARNING: No matching mzML for {filter_idxml}")
            continue

        psms_out = str(library_dir / f"{original_stem}_psms.tsv")
        peaks_out = str(library_dir / f"{original_stem}_peaks.tsv")

        convert_cmd = [
            "easypqp", "convert",
            "--pepxml", filter_idxml,
            "--spectra", matching_mzml,
            "--psms", psms_out,
            "--peaks", peaks_out,
        ]

        with st.spinner(f"Converting {Path(filter_idxml).name}..."):
            if not self.executor.run_command(convert_cmd):
                self.logger.log(f"WARNING: EasyPQP convert failed for {filter_idxml}")
                continue

        psms_files.append(psms_out)
        peaks_files.append(peaks_out)

    if psms_files:
        # Step 2: Build merged library
        library_pqp = str(library_dir / "spectral_library.pqp")

        library_cmd = ["easypqp", "library", "--out", library_pqp]

        if not self.params.get("library-use-fdr", False):
            library_cmd.append("--nofdr")
        else:
            fdr = self.params.get("library-psm-fdr", 0.01)
            library_cmd.extend(["--psm_fdr_threshold", str(fdr)])

        for psms, peaks in zip(psms_files, peaks_files):
            library_cmd.extend([psms, peaks])

        with st.spinner("Building spectral library..."):
            if not self.executor.run_command(library_cmd):
                self.logger.log("ERROR: EasyPQP library failed")
                st.error("Library generation failed.")
            else:
                final_library = library_pqp

                # Step 3: Generate decoys (optional)
                if self.params.get("library-generate-decoys", True):
                    method = self.params.get("library-decoy-method", "shuffle")
                    library_decoy = str(library_dir / "spectral_library_decoy.pqp")

                    decoy_cmd = [
                        "easypqp", "openswath_decoy_generator",
                        "--in", library_pqp,
                        "--out", library_decoy,
                        "--method", method,
                    ]

                    with st.spinner("Generating decoys..."):
                        if self.executor.run_command(decoy_cmd):
                            final_library = library_decoy
                        else:
                            self.logger.log("WARNING: Decoy generation failed")

                # Step 4: Convert to TSV
                library_tsv = str(library_dir / "spectral_library.tsv")

                tsv_cmd = [
                    "TargetedFileConverter",
                    "-in", final_library,
                    "-out", library_tsv,
                ]

                with st.spinner("Converting to TSV..."):
                    if self.executor.run_command(tsv_cmd):
                        self.logger.log(f"Library created: {library_tsv}")
                        st.success("Spectral library created successfully!")
                    else:
                        self.logger.log("WARNING: TSV conversion failed")
    else:
        self.logger.log("ERROR: No files converted for library")
        st.error("Library generation failed: No files converted.")
```

---

## 3. Results Page Changes (`content/results_filtered.py`)

### Insert Before Line 75 (`st.markdown("---")`)

```python
# Spectral Library Download Section
st.markdown("---")
st.subheader("Spectral Library")

library_dir = workflow_dir / "results" / "library"
library_tsv = library_dir / "spectral_library.tsv"

if library_tsv.exists():
    st.success("Spectral library available for download.")
    with open(library_tsv, "rb") as f:
        st.download_button(
            label="Download Spectral Library (TSV)",
            data=f,
            file_name="spectral_library.tsv",
            mime="text/tab-separated-values",
            use_container_width=True,
        )
elif library_dir.exists():
    st.warning("Library generation was enabled but TSV file not found.")
else:
    st.info("Spectral library generation was not enabled for this workflow run.")
```

---

## 4. Default Parameters (`default-parameters.json`)

Add these entries:

```json
{
    "generate-library": false,
    "library-use-fdr": false,
    "library-psm-fdr": 0.01,
    "library-generate-decoys": true,
    "library-decoy-method": "shuffle"
}
```

---

## Directory Structure (New)

```
workflow_dir/results/
├── filter_results/
│   └── {sample}_filter.idXML
├── library/                      # NEW
│   ├── {sample}_psms.tsv
│   ├── {sample}_peaks.tsv
│   ├── spectral_library.pqp
│   ├── spectral_library_decoy.pqp
│   └── spectral_library.tsv      # Final output
└── quant_results/
```

---

## Error Handling

| Scenario | Behavior |
|----------|----------|
| EasyPQP not installed | `run_command()` fails, logs warning, workflow continues |
| No matching mzML | Skip that file, log warning, continue with others |
| Convert fails for a file | Skip that file, continue with others |
| No files converted | Show error, skip library step |
| Decoy generation fails | Use library without decoys, warn user |
| TSV conversion fails | Log warning, PQP still available |

---

## Verification

1. Run app: `streamlit run app.py`
2. Configure workflow with library generation enabled
3. Run workflow
4. Check `results/library/` directory for output files
5. Navigate to Filtered PSMs results page
6. Verify download button appears and downloads valid TSV
