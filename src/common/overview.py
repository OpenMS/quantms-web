"""Run overview: key figures and a per-file table of a finished workflow run.

Reads only what the workflow already wrote: the input mzML files (MS2 spectra
are counted without loading peaks), the Percolator and IDFilter idXMLs, and the
run-time log. It makes no Streamlit calls, so the page caches it and tests call it.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pyopenms import AASequence, MSExperiment, MzMLFile, PeakFileOptions, SpectrumSettings

from src.common.contaminants import CONTAMINANT_PREFIX
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

# Columns the identification diagnosis (src/common/id_diagnosis.py) reads besides
# PER_FILE_COLUMNS: the share of MS2 spectra stored as profile data, the share of
# target PSMs on a contaminant (trypsin left out, it is in every digest) and the
# share whose precursor error lies at the edge of the searched window
DIAGNOSIS_COLUMNS = ["Profile MS2 (%)", "Contaminant PSMs (%)", "At window edge (%)"]

# A PSM is at the window edge when its error is beyond this share of the tolerance
EDGE_OF_WINDOW = 0.9


def ms2_summary(mzml_path) -> tuple[int | None, float | None]:
    """MS2 spectra of an mzML and the share (%) stored as profile data, read without peaks.

    (None, None) if the file is missing or unreadable; the share is None
    without MS2 spectra.
    """
    try:
        options = PeakFileOptions()
        options.setMSLevels([2])
        options.setFillData(False)
        mzml = MzMLFile()
        mzml.setOptions(options)
        exp = MSExperiment()
        mzml.load(str(mzml_path), exp)
    except Exception:  # missing or unreadable file: shown as empty, not an error
        return None, None
    n = exp.getNrSpectra()
    if not n:
        return 0, None
    # Spectra of unknown type count as centroided
    profile_type = SpectrumSettings.SpectrumType.PROFILE
    profile = sum(1 for spec in exp.getSpectra() if spec.getType() == profile_type)
    return n, round(100 * profile / n, 1)


def count_ms2_spectra(mzml_path) -> int | None:
    """Number of MS2 spectra in an mzML, read without peak data; None if unreadable."""
    return ms2_summary(mzml_path)[0]


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


def _is_contaminant(hit) -> bool:
    """Whether a hit's proteins are all contaminants other than trypsin."""
    accessions = [e.getProteinAccession() for e in hit.getPeptideEvidences()]
    accessions = [a.decode() if isinstance(a, bytes) else str(a) for a in accessions]
    return bool(accessions) and all(
        a.startswith(CONTAMINANT_PREFIX) and "TRYP_" not in a for a in accessions
    )


PSM_COLUMNS = ["sequence", "decoy", "ppm", "contaminant"]


def read_psms(idxml_path) -> pd.DataFrame:
    """Top hit of every spectrum in an idXML: sequence, decoy flag, precursor error, contaminant flag."""
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
            "contaminant": _is_contaminant(hit),
        })
    return pd.DataFrame(rows, columns=PSM_COLUMNS)


def search_times(log_path, mzml_files) -> dict:
    """Seconds the database search (ProSE) took per input file, from the latest run in the run-time log."""
    path = Path(log_path)
    if not path.is_file():
        return {}
    text = path.read_text(encoding="utf-8", errors="ignore")
    times = {}
    blocks = re.findall(
        r"Process finished:\n((?:ProSE|CometAdapter)\b.*?)\nTotal time to run command: ([\d.]+) seconds",
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


def _share(mask: pd.Series) -> float | None:
    """Percentage of True in a boolean Series, None when it is empty."""
    return round(100 * float(mask.astype(bool).mean()), 1) if len(mask) else None


def run_overview(workflow_dir, mzml_files, precursor_ppm: float | None = None) -> tuple[pd.DataFrame, dict]:
    """Per-file table and run-wide key figures of a finished run.

    The table has PER_FILE_COLUMNS, then DIAGNOSIS_COLUMNS. PSM counts are at
    the PSM FDR set for the run (the IDFilter output). Files the run did not
    reach show empty cells. ``precursor_ppm`` is the searched precursor
    tolerance in ppm (the wider side), for the share of PSMs at its edge;
    without it (a tolerance in Da) that column stays empty.
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
        ms2, profile = ms2_summary(mzml)
        psms = read_psms(filtered[name]) if name in filtered else None
        # astype(bool): an empty file gives an object column, and indexing with an
        # empty object Series selects columns instead of rows
        targets = psms[~psms["decoy"].astype(bool)] if psms is not None else None
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
            "Profile MS2 (%)": profile,
            "Contaminant PSMs (%)": _share(targets["contaminant"]) if targets is not None else None,
            "At window edge (%)": (
                _share(targets["ppm"].dropna().abs() >= EDGE_OF_WINDOW * precursor_ppm)
                if targets is not None and precursor_ppm else None
            ),
        })
    per_file = pd.DataFrame(rows, columns=PER_FILE_COLUMNS + DIAGNOSIS_COLUMNS)

    targets = pd.concat(all_psms) if all_psms else pd.DataFrame(columns=PSM_COLUMNS)
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
