"""Run overview: key figures and a per-file table of a finished workflow run.

Reads only what the workflow already wrote: the input mzML files (MS2 spectra
are counted without loading peaks), the Percolator and IDFilter idXMLs, and the
run-time log. It makes no Streamlit calls, so the page caches it and tests call it.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pyopenms import AASequence, MSExperiment, MzMLFile, PeakFileOptions

from src.common.results_helpers import extract_filename_from_idxml, load_idxml

# Mass difference between the 13C and 12C isotopes; precursor picks up to
# three isotope peaks off the monoisotope are corrected before the error is taken
C13_DELTA = 1.0033548378
MAX_ISOTOPE_ERROR = 3

PER_FILE_COLUMNS = [
    "File",
    "MS2 spectra",
    "Target PSMs",
    "Decoy PSMs",
    "ID rate (%)",
    "Peptides",
    "Mass error (ppm)",
    "Percolator",
    "Search (s)",
]


def count_ms2_spectra(mzml_path) -> int | None:
    """Number of MS2 spectra in an mzML, read without peak data; None if unreadable."""
    try:
        options = PeakFileOptions()
        options.setMSLevels([2])
        options.setFillData(False)
        mzml = MzMLFile()
        mzml.setOptions(options)
        exp = MSExperiment()
        mzml.load(str(mzml_path), exp)
        return exp.getNrSpectra()
    except Exception:  # missing or unreadable file: shown as empty, not an error
        return None


def precursor_error_ppm(sequence: AASequence, charge: int, observed_mz: float) -> float | None:
    """Precursor m/z error in ppm, after correcting a 13C isotope pick of up to three peaks."""
    if charge <= 0:
        return None
    theoretical = sequence.getMZ(charge)
    errors = [
        (observed_mz - k * C13_DELTA / charge - theoretical) / theoretical * 1e6
        for k in range(MAX_ISOTOPE_ERROR + 1)
    ]
    return min(errors, key=abs)


def read_psms(idxml_path) -> pd.DataFrame:
    """Top hit of every spectrum in an idXML: sequence, decoy flag and precursor error."""
    _, peptides = load_idxml(idxml_path)
    rows = []
    for pep in peptides:
        hits = pep.getHits()
        if not hits:
            continue
        hit = hits[0]
        target_decoy = hit.getMetaValue("target_decoy") if hit.metaValueExists("target_decoy") else ""
        if isinstance(target_decoy, bytes):
            target_decoy = target_decoy.decode()
        sequence = hit.getSequence()
        rows.append({
            "sequence": sequence.toString(),
            "decoy": str(target_decoy).startswith("decoy"),
            "ppm": precursor_error_ppm(sequence, hit.getCharge(), pep.getMZ()),
        })
    return pd.DataFrame(rows, columns=["sequence", "decoy", "ppm"])


def search_times(log_path, mzml_files) -> dict:
    """Seconds CometAdapter took per input file, from the latest run in the run-time log."""
    path = Path(log_path)
    if not path.is_file():
        return {}
    text = path.read_text(encoding="utf-8", errors="ignore")
    times = {}
    blocks = re.findall(
        r"Process finished:\n(CometAdapter\b.*?)\nTotal time to run command: ([\d.]+) seconds",
        text,
        flags=re.DOTALL,
    )
    for command, seconds in blocks:
        match = re.search(r"\s-in\s+(\S+)", command)
        if match:
            # Later runs append to the log, so the last block per file wins
            times[Path(match.group(1)).name] = float(seconds)
    return {Path(f).name: times[Path(f).name] for f in mzml_files if Path(f).name in times}


def _by_input(directory: Path, names: list) -> dict:
    """idXMLs in ``directory`` keyed by the input mzML file name they came from."""
    if not directory.is_dir():
        return {}
    return {extract_filename_from_idxml(p, names): p for p in sorted(directory.glob("*.idXML"))}


def run_overview(workflow_dir, mzml_files) -> tuple[pd.DataFrame, dict]:
    """Per-file table (PER_FILE_COLUMNS) and run-wide key figures of a finished run.

    PSM counts are at the PSM FDR set for the run (the IDFilter output). Files
    the run did not reach show empty cells.
    """
    workflow_dir = Path(workflow_dir)
    results = workflow_dir / "results"
    names = [Path(f).name for f in mzml_files]
    filtered = _by_input(results / "psm_filter", names)
    rescored = _by_input(results / "percolator_results", names)
    times = search_times(workflow_dir / "logs" / "commands-and-run-times.log", mzml_files)

    rows = []
    all_psms = []
    for mzml in mzml_files:
        name = Path(mzml).name
        ms2 = count_ms2_spectra(mzml)
        psms = read_psms(filtered[name]) if name in filtered else None
        targets = psms[~psms["decoy"]] if psms is not None else None
        if targets is not None:
            all_psms.append(targets)
        rows.append({
            "File": name,
            "MS2 spectra": ms2,
            "Target PSMs": len(targets) if targets is not None else None,
            "Decoy PSMs": int(psms["decoy"].sum()) if psms is not None else None,
            "ID rate (%)": round(100 * len(targets) / ms2, 1) if targets is not None and ms2 else None,
            "Peptides": targets["sequence"].nunique() if targets is not None else None,
            "Mass error (ppm)": (
                round(float(targets["ppm"].median()), 2)
                if targets is not None and targets["ppm"].notna().any() else None
            ),
            "Percolator": name in rescored,
            "Search (s)": round(times[name]) if name in times else None,
        })
    per_file = pd.DataFrame(rows, columns=PER_FILE_COLUMNS)

    targets = pd.concat(all_psms) if all_psms else pd.DataFrame(columns=["sequence", "decoy", "ppm"])
    ms2_total = per_file["MS2 spectra"].dropna().sum() if per_file["MS2 spectra"].notna().any() else None
    key_figures = {
        "Target PSMs": len(targets),
        "Peptides": targets["sequence"].nunique(),
        "PSMs per MS2 (%)": round(100 * len(targets) / ms2_total, 1) if ms2_total else None,
        "Median precursor error (ppm)": (
            round(float(np.nanmedian(targets["ppm"].astype(float))), 2) if targets["ppm"].notna().any() else None
        ),
        "Files": len(mzml_files),
    }
    return per_file, key_figures
