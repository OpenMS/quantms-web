"""Why a file has few identifications: the files of a run that fall short, and likely causes.

Reads the per-file table of src/common/overview.py (PER_FILE_COLUMNS and
DIAGNOSIS_COLUMNS) and the run's settings. Streamlit-free, so tests call it.

A file falls short when its ID rate (target PSMs at the PSM FDR per MS2
spectrum) is below what its kind of search usually reaches, or below a third
of the median of the run's files. For each such file it names up to three
causes, the most likely first.
"""
import pandas as pd

from src.common.results_helpers import _tool_params, search_tolerances

# ID rates (%) a search of this kind usually reaches on good data
EXPECTED_ID_RATE = {
    "tryptic": 15.0,
    "ion trap": 8.0,
    "nonspecific": 3.0,
}
# A file below this share of the run's median ID rate falls short
BATCH_SHARE = 1 / 3
# The run's median counts once it has this many files with an ID rate
MIN_FILES_FOR_MEDIAN = 3
# Shares (%) above which a cause is named
PROFILE_SHARE = 10.0
CONTAMINANT_SHARE = 30.0
EDGE_SHARE = 10.0
# A median precursor error beyond this share of the tolerance is an offset
OFFSET_SHARE = 0.5
SMALL_RUN_PSMS = 300
MAX_CAUSES = 3

DIAGNOSIS_COLUMNS = ["File", "ID rate (%)", "Expected (%)", "Likely causes"]


def search_kind(workflow_params: dict) -> str:
    """'nonspecific' (no enzyme), 'ion trap' (fragment tolerance in Da) or 'tryptic'."""
    tmt = workflow_params.get("analysis-mode", "LFQ") == "TMT"
    prose = _tool_params(workflow_params, "ProSE-TMT" if tmt else "ProSE")
    enzyme = str(prose.get("Search:enzyme", "Trypsin/P")).lower()
    if enzyme in ("unspecific cleavage", "no cleavage", "unspecific", "none"):
        return "nonspecific"
    if not search_tolerances(workflow_params)["fragment_ppm"]:
        return "ion trap"
    return "tryptic"


def _value(row, column):
    """A cell of the per-file table, None when it is empty."""
    value = row.get(column)
    return None if value is None or value != value else value


def _causes(row, tolerances: dict) -> list:
    """Likely causes of one file's low ID rate, the most likely first."""
    causes = []  # (weight, text)
    ms2 = _value(row, "MS2 spectra")
    if not ms2:
        causes.append((1.0, "The file has no MS2 spectra to search."))
    profile = _value(row, "Profile MS2 (%)")
    if profile is not None and profile > PROFILE_SHARE:
        causes.append((0.9, f"{profile:g} % of its MS2 spectra are profile data; the search expects "
                            "centroided spectra (peak picking when converting the raw file)."))
    error = _value(row, "Mass error (ppm)")
    if tolerances["precursor_unit"] == "ppm" and tolerances["precursor"] < 100:
        window = tolerances["precursor"]
        if error is not None and abs(error) > OFFSET_SHARE * window:
            causes.append((0.7, f"Its median precursor error is {error:+g} ppm in a ±{window:g} ppm "
                                "window: the instrument's calibration is off, so a wider precursor "
                                "tolerance would find more."))
        edge = _value(row, "At window edge (%)")
        offset = error is not None and abs(error) > OFFSET_SHARE * window
        # An offset puts PSMs at the edge too; naming it once is enough
        if edge is not None and edge > EDGE_SHARE and not offset:
            causes.append((0.7, f"{edge:g} % of its PSMs lie at the edge of the ±{window:g} ppm window: "
                                "the precursor tolerance is likely too narrow."))
    contaminants = _value(row, "Contaminant PSMs (%)")
    if contaminants is not None and contaminants > CONTAMINANT_SHARE:
        causes.append((0.6, f"{contaminants:g} % of its PSMs are on contaminants (keratins, BSA, ...): "
                            "the sample itself gave little."))
    if _value(row, "Percolator") is False and ms2:
        causes.append((0.5, "Percolator did not rescore it, so the search scores alone set the FDR."))
    targets = _value(row, "Target PSMs")
    if targets is not None and targets < SMALL_RUN_PSMS and ms2:
        causes.append((0.3, f"With {targets:,} target PSMs the run is small, and the FDR estimate "
                            "removes a larger share of a small run."))
    causes.sort(key=lambda c: -c[0])
    return [text for _, text in causes[:MAX_CAUSES]]


def diagnose_identifications(per_file: pd.DataFrame, workflow_params: dict, all_files: bool = False) -> pd.DataFrame:
    """Files of a run whose ID rate falls short, with the rate expected and the likely causes.

    One row per such file (DIAGNOSIS_COLUMNS; the causes joined by spaces),
    empty when every file reaches the expected rate. With ``all_files`` the
    other files follow, their causes empty. A file the run did not search (no
    ID rate and no MS2 count) is left out.
    """
    expected = EXPECTED_ID_RATE[search_kind(workflow_params)]
    tolerances = search_tolerances(workflow_params)
    rates = pd.to_numeric(per_file["ID rate (%)"], errors="coerce")
    median = rates.median() if rates.notna().sum() >= MIN_FILES_FOR_MEDIAN else None
    rows, fine = [], []
    for (_, row), rate in zip(per_file.iterrows(), rates):
        no_ms2 = _value(row, "MS2 spectra") == 0
        if rate != rate and not no_ms2:
            continue
        low = no_ms2 or rate < expected
        outlier = median is not None and not no_ms2 and rate < BATCH_SHARE * median
        if not (low or outlier):
            fine.append({"File": row["File"], "ID rate (%)": float(rate), "Expected (%)": expected,
                         "Likely causes": ""})
            continue
        causes = _causes(row, tolerances)
        if outlier and not low:
            causes.append(f"Well below the other files of the run (median {median:g} %).")
        if not causes:
            causes = ["None of the checks explains it: check that the database, enzyme and "
                      "modifications fit this sample."]
        rows.append({
            "File": row["File"],
            "ID rate (%)": None if rate != rate else float(rate),
            "Expected (%)": expected,
            "Likely causes": " ".join(causes[:MAX_CAUSES]),
        })
    return pd.DataFrame(rows + (fine if all_files else []), columns=DIAGNOSIS_COLUMNS)
