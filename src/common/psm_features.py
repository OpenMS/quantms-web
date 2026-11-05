"""Every score and feature of a PSM, read from its idXML, with readable names.

ProSE writes its score and the features of each match under plain names (runs
made before the switch to ProSE carry Comet's scores, under a PSI-MS accession
and a ``COMET:`` name), and Percolator adds its score, q-value and PEP as PSI-MS
accessions. This module
names them, drops the duplicates and keeps any value it does not know under
its own key, so a new tool version never hides a feature.
Streamlit-free, so it can be tested on its own.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.common.results_helpers import load_idxml

# key in the idXML -> (group, label, help)
KNOWN_FEATURES: dict[str, tuple[str, str, str]] = {
    "MS:1001492": ("Percolator", "Percolator score", "Percolator's SVM score: higher is more confident."),
    "MS:1001491": ("Percolator", "q-value", "Lowest FDR at which this PSM is accepted."),
    "MS:1001493": ("Percolator", "Posterior error probability (PEP)",
                   "Probability that this single PSM is wrong."),
    "ln(hyperscore)": ("ProSE", "ln(hyperscore)",
                       "ProSE's score: the matched b and y ions, their number and intensity. Higher is better."),
    "delta_score": ("ProSE", "Δ score", "Score gap to the next-best peptide for this spectrum."),
    "hyperscore_zscore": ("ProSE", "Score z-score", "The score against those of all candidate peptides of the spectrum."),
    "ln_num_candidates": ("ProSE", "ln(candidate peptides)", "Natural log of the peptides scored for this spectrum."),
    "num_matched_peaks": ("ProSE", "Matched peaks", "Peaks explained by a fragment ion of the peptide."),
    "matched_prefix_ions": ("ProSE", "Matched b ions", "Fragment ions from the peptide's start found in the spectrum."),
    "matched_suffix_ions": ("ProSE", "Matched y ions", "Fragment ions from the peptide's end found in the spectrum."),
    "matched_prefix_ions_fraction": ("ProSE", "Fraction of b ions matched", "Matched / possible b ions."),
    "matched_suffix_ions_fraction": ("ProSE", "Fraction of y ions matched", "Matched / possible y ions."),
    "longest_peptide_ion_sequence": ("ProSE", "Longest ion ladder", "Most consecutive matched b or y ions."),
    "complementary_ions_fraction": ("ProSE", "Complementary ion pairs", "Share of cleavage sites seen as both a b and a y ion."),
    "matched_ion_current": ("ProSE", "Matched intensity", "Intensity of the matched peaks."),
    "matched_ion_current_fraction": ("ProSE", "Fraction of intensity matched", "Matched / total peak intensity."),
    "fragment_mz_error_median_ppm": ("ProSE", "Median fragment error (ppm)", "Median m/z error of the matched fragment ions."),
    "precursor_mz_error_ppm": ("ProSE", "Precursor error (ppm)", "m/z error of the precursor, after the isotope error."),
    "MS:1002252": ("Comet", "XCorr", "Cross-correlation of the measured and theoretical spectrum."),
    "MS:1002253": ("Comet", "ΔCn", "Relative XCorr gap to the next-best peptide for this spectrum."),
    "MS:1002255": ("Comet", "Sp score", "Comet's preliminary score."),
    "MS:1002256": ("Comet", "Sp rank", "Rank of this peptide by Sp score."),
    "MS:1002257": ("Comet", "E-value", "Expected number of random matches scoring this well: lower is better."),
    "MS:1002258": ("Comet", "Matched fragment ions", "Theoretical fragment ions found in the spectrum."),
    "MS:1002259": ("Comet", "Total fragment ions", "Theoretical fragment ions of the peptide."),
    "COMET:deltaLCn": ("Comet", "ΔLCn", "XCorr gap to the last-ranked peptide for this spectrum."),
    "COMET:lnExpect": ("Comet", "ln(E-value)", "Natural log of the E-value, as Percolator reads it."),
    "COMET:lnRankSP": ("Comet", "ln(Sp rank)", "Natural log of the Sp rank."),
    "COMET:lnNumSP": ("Comet", "ln(peptides scored)", "Natural log of the number of peptides Comet scored for this spectrum."),
    "COMET:IonFrac": ("Comet", "Fraction of ions matched", "Matched / total fragment ions."),
    "num_matched_peptides": ("Search", "Candidate peptides", "Peptides within the precursor tolerance."),
    "isotope_error": ("Search", "Isotope error", "Precursor picked off the monoisotopic peak by this many 13C."),
    "target_decoy": ("Search", "Target or decoy", "Target: from the real database. Decoy: from the reversed one."),
}

# Same value as a known key, written again under another name
DUPLICATES = {
    "percolator_score", "percolator_q_value", "percolator_pep",
    "COMET:xcorr", "COMET:deltaCn", "COMET:spscore", "COMET:sprank", "Comet:lnrSp", "Comet:IonFrac", "expect",
}
# Not features
SKIPPED = {"protein_references", "spectrum_reference", "fragment_annotation", "scan_index"}

GROUP_ORDER = ["Percolator", "ProSE", "Comet", "Search", "Other"]


def _plain(value):
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, float):
        return float(f"{value:.6g}")
    return value


def hit_features(hit) -> list[dict]:
    """Rows (group, feature, value, help) of one PeptideHit, Percolator first.

    The hit's main score is not listed on its own: it is always one of these
    features too (ProSE's ln(hyperscore), or Percolator's PEP after rescoring).
    """
    known = list(KNOWN_FEATURES)
    rows = []
    keys: list = []
    hit.getKeys(keys)
    for raw in keys:
        key = raw.decode() if isinstance(raw, bytes) else str(raw)
        if key in DUPLICATES or key in SKIPPED:
            continue
        group, label, help_text = KNOWN_FEATURES.get(key, ("Other", key, ""))
        rows.append({"group": group, "feature": label, "value": _plain(hit.getMetaValue(raw)), "help": help_text,
                     "_rank": known.index(key) if key in KNOWN_FEATURES else len(known)})
    rows.sort(key=lambda r: (GROUP_ORDER.index(r["group"]), r["_rank"]))
    for row in rows:
        del row["_rank"]
    return rows


def features_by_id(idxml: Path) -> dict[int, list[dict]]:
    """Feature rows of every hit in an idXML, keyed by the ``id_idx`` that ``parse_idxml`` gives it."""
    _, peptides = load_idxml(idxml)
    features: dict[int, list[dict]] = {}
    for pep in peptides:
        for hit in pep.getHits():
            features[len(features)] = hit_features(hit)
    return features


def features_frame(rows: list[dict]) -> pd.DataFrame:
    """Rows as a table for display, values as text so mixed types line up."""
    frame = pd.DataFrame(rows, columns=["group", "feature", "value", "help"])
    frame["value"] = frame["value"].map(lambda v: f"{v:g}" if isinstance(v, float) else str(v))
    return frame.rename(columns={"group": "Source", "feature": "Feature", "value": "Value", "help": "Meaning"})
