"""AlphaPeptDeep prediction features for Percolator.

For every PSM of a search result this adds three features that a search
engine score cannot see:

- ``peptdeep:spectral_angle``: normalised spectral contrast angle between the
  measured spectrum and AlphaPeptDeep's predicted b/y intensities (1 identical,
  0 nothing in common);
- ``peptdeep:cosine``: cosine similarity of the same two vectors;
- ``peptdeep:delta_rt``: absolute difference in seconds between the measured
  retention time and the predicted one, after a linear calibration of the
  predictions on the file's confident target PSMs.

The collision energy is calibrated per file: the NCE whose predictions best
match a sample of confident target PSMs (highest median spectral angle) is
used for all PSMs of that file. Targets and decoys get features the same way,
so Percolator can learn from them. The feature names are added to the
search's ``extra_features`` list, which PercolatorAdapter reads.

Streamlit-free; ``src/python-tools/peptdeep_rescore.py`` runs it as its own
process in the workflow, so the model's memory is freed afterwards.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, List, Optional

import numpy as np
import pandas as pd

from src.common.ms2_prediction import (
    DEFAULT_INSTRUMENT,
    FRAG_TYPES,
    instrument_of,
    match_peaks,
    spectral_angle,
    to_peptdeep,
)

FEATURES = ["peptdeep:spectral_angle", "peptdeep:cosine", "peptdeep:delta_rt"]
NCE_GRID = [20.0, 23.0, 26.0, 28.0, 30.0, 32.0, 35.0, 38.0]
CALIBRATION_PSMS = 300
CONFIDENT_EXPECT = 0.01  # e-value of the PSMs used to calibrate NCE and RT, where the search reports one
ALLOWED = set("ACDEFGHIKLMNPQRSTVWY")


def cosine(observed: np.ndarray, predicted: np.ndarray) -> float:
    a, b = np.asarray(observed, float), np.asarray(predicted, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return float(np.clip(a @ b / (na * nb), 0.0, 1.0))


def psm_table(peptides) -> pd.DataFrame:
    """One row per PSM (top hit of each spectrum) with what prediction and scoring need."""
    rows = []
    for i, pep in enumerate(peptides):
        hits = pep.getHits()
        if not hits:
            continue
        hit = hits[0]
        peptide = hit.getSequence().toString()
        converted = to_peptdeep(peptide)
        expect = hit.getMetaValue("expect") if hit.metaValueExists("expect") else None
        # Best first: the e-value where there is one, else the search score (ProSE: ln(hyperscore))
        score = float(hit.getScore())
        rank = float(expect) if expect is not None else (-score if pep.isHigherScoreBetter() else score)
        td = hit.getMetaValue("target_decoy") if hit.metaValueExists("target_decoy") else ""
        td = td.decode() if isinstance(td, bytes) else str(td)
        ref = pep.getMetaValue("spectrum_reference") if pep.metaValueExists("spectrum_reference") else ""
        rows.append({
            "pep_index": i,
            "spectrum_reference": ref.decode() if isinstance(ref, bytes) else str(ref),
            "rt": float(pep.getRT()),
            "charge": int(hit.getCharge()),
            "decoy": td.startswith("decoy"),
            "expect": float(expect) if expect is not None else np.nan,
            "rank": rank,
            "sequence": converted[0] if converted else None,
            "mods": converted[1] if converted else None,
            "mod_sites": converted[2] if converted else None,
        })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    ok = df["sequence"].notna()
    ok &= df["sequence"].fillna("").map(lambda s: 3 <= len(s) <= 60 and not set(s) - ALLOWED)
    ok &= df["charge"].between(1, 6)
    df["predictable"] = ok
    return df


def predict(model, psms: pd.DataFrame, nce: float, instrument: str, items=("ms2", "rt")) -> dict:
    """peptdeep predictions for the predictable PSMs (rt_pred and fragment tables)."""
    precursors = psms.loc[psms["predictable"], ["sequence", "mods", "mod_sites", "charge"]].copy()
    precursors["nce"] = float(nce)
    precursors["instrument"] = instrument
    precursors["_row"] = precursors.index
    precursors = precursors.reset_index(drop=True)
    result = model.predict_all(precursors, predict_items=list(items), frag_types=FRAG_TYPES, multiprocessing=False)
    return result


def fragment_vectors(result: dict, row: int, charge: int):
    """(m/z, intensity) arrays of one precursor's predicted ions, fragment charges below its charge."""
    precursors = result["precursor_df"]
    start, stop = int(precursors["frag_start_idx"].iat[row]), int(precursors["frag_stop_idx"].iat[row])
    mzs, intensities = [], []
    for frag_type in FRAG_TYPES:
        z = int(frag_type.split("_z")[1])
        if z > 1 and z >= charge:
            continue
        mzs.append(result["fragment_mz_df"][frag_type].to_numpy()[start:stop])
        intensities.append(result["fragment_intensity_df"][frag_type].to_numpy()[start:stop])
    mz, inten = np.concatenate(mzs), np.concatenate(intensities)
    keep = (mz > 0) & (inten > 0)
    return mz[keep], inten[keep]


def similarity(result: dict, psms: pd.DataFrame, spectra: Dict[str, tuple], tolerance: float, ppm: bool) -> pd.DataFrame:
    """Spectral angle and cosine of every predicted PSM (indexed like ``psms``)."""
    precursors = result["precursor_df"]
    out = {}
    for row in range(len(precursors)):
        idx = int(precursors["_row"].iat[row])
        psm = psms.loc[idx]
        spectrum = spectra.get(psm["spectrum_reference"])
        if spectrum is None:
            continue
        mz, inten = fragment_vectors(result, row, int(psm["charge"]))
        if len(mz) == 0:
            continue
        observed = match_peaks(spectrum[0], np.sqrt(spectrum[1]), mz, tolerance, ppm)
        predicted = np.sqrt(inten)
        out[idx] = (spectral_angle(observed, predicted), cosine(observed, predicted))
    return pd.DataFrame.from_dict(out, orient="index", columns=["spectral_angle", "cosine"])


def calibration_set(psms: pd.DataFrame, n: int = CALIBRATION_PSMS) -> pd.DataFrame:
    """Confident target PSMs (best e-values, or best search scores) to calibrate NCE and RT on."""
    order = "rank" if "rank" in psms else "expect"
    targets = psms[psms["predictable"] & ~psms["decoy"]].sort_values(order)
    confident = targets[targets["expect"] <= CONFIDENT_EXPECT]
    if len(confident) < 30:
        confident = targets.head(max(30, n // 3))
    return confident.head(n)


def calibrate_nce(model, psms: pd.DataFrame, spectra, instrument: str, tolerance: float, ppm: bool,
                  grid: List[float] = NCE_GRID) -> float:
    """The NCE whose predictions match the confident PSMs best (highest median spectral angle)."""
    sample = calibration_set(psms)
    if sample.empty:
        return 30.0
    best, best_score = 30.0, -1.0
    for nce in grid:
        result = predict(model, sample, nce, instrument, items=("ms2",))
        scores = similarity(result, sample, spectra, tolerance, ppm)
        score = float(scores["spectral_angle"].median()) if not scores.empty else -1.0
        if score > best_score:
            best, best_score = nce, score
    return best


def calibrate_rt(predicted: pd.Series, observed: pd.Series, confident: pd.Index) -> Callable:
    """Linear map from peptdeep's normalised RT to seconds, fitted on the confident PSMs.

    The fit is repeated once without the worst 20 % of residuals, so a few
    wrong identifications do not tilt it.
    """
    x = predicted.loc[confident].to_numpy(float)
    y = observed.loc[confident].to_numpy(float)
    if len(x) < 10 or np.ptp(x) == 0:
        span = max(float(observed.max() - observed.min()), 1.0)
        return lambda p: observed.min() + np.asarray(p, float) * span
    slope, intercept = np.polyfit(x, y, 1)
    residual = np.abs(y - (slope * x + intercept))
    keep = residual <= np.quantile(residual, 0.8)
    if keep.sum() >= 10:
        slope, intercept = np.polyfit(x[keep], y[keep], 1)
    return lambda p: slope * np.asarray(p, float) + intercept


def features(model, psms: pd.DataFrame, spectra: Dict[str, tuple], instrument: str, nce: float,
             tolerance: float, ppm: bool) -> pd.DataFrame:
    """The three features of every PSM; unpredictable PSMs get the worst values of the file."""
    result = predict(model, psms, nce, instrument)
    sim = similarity(result, psms, spectra, tolerance, ppm)
    precursors = result["precursor_df"]
    rt_pred = pd.Series(precursors["rt_pred"].to_numpy(float), index=precursors["_row"].to_numpy())
    observed = psms["rt"]
    confident = calibration_set(psms).index.intersection(rt_pred.index)
    to_seconds = calibrate_rt(rt_pred, observed.loc[rt_pred.index], confident)
    delta = (observed.loc[rt_pred.index] - pd.Series(to_seconds(rt_pred.to_numpy()), index=rt_pred.index)).abs()
    out = pd.DataFrame(index=psms.index)
    out["peptdeep:spectral_angle"] = sim["spectral_angle"]
    out["peptdeep:cosine"] = sim["cosine"]
    out["peptdeep:delta_rt"] = delta
    worst_rt = float(delta.max()) if len(delta) else float(observed.max() - observed.min())
    return out.fillna({"peptdeep:spectral_angle": 0.0, "peptdeep:cosine": 0.0, "peptdeep:delta_rt": worst_rt})


def read_spectra(mzml_path) -> Dict[str, tuple]:
    """Native ID -> (m/z, intensity) of every MS2 spectrum."""
    from pyopenms import MSExperiment, MzMLFile, PeakFileOptions

    options = PeakFileOptions()
    options.setMSLevels([2])
    mzml = MzMLFile()
    mzml.setOptions(options)
    exp = MSExperiment()
    mzml.load(str(mzml_path), exp)
    spectra = {}
    for spectrum in exp:
        mz, inten = spectrum.get_peaks()
        spectra[spectrum.getNativeID()] = (np.asarray(mz, float), np.asarray(inten, float))
    return spectra


def add_features(idxml_in, mzml, idxml_out, model, tolerance: float, ppm: bool,
                 instrument: Optional[str] = None, log=print) -> dict:
    """Annotate an idXML's PSMs with the features and register them for Percolator."""
    from pyopenms import IdXMLFile
    from src.common.results_helpers import load_idxml

    proteins, peptides = load_idxml(idxml_in)
    psms = psm_table(peptides)
    spectra = read_spectra(mzml)
    instrument = instrument or instrument_of(mzml) or DEFAULT_INSTRUMENT
    nce = calibrate_nce(model, psms, spectra, instrument, tolerance, ppm) if not psms.empty else 30.0
    values = features(model, psms, spectra, instrument, nce, tolerance, ppm) if not psms.empty else pd.DataFrame()

    for idx, row in values.iterrows():
        pep = peptides[int(psms.at[idx, "pep_index"])]
        hits = list(pep.getHits())
        for name in FEATURES:
            hits[0].setMetaValue(name, float(row[name]))
        pep.setHits(hits)
    for prot in proteins:
        search = prot.getSearchParameters()
        extra = search.getMetaValue("extra_features") if search.metaValueExists("extra_features") else ""
        extra = extra.decode() if isinstance(extra, bytes) else str(extra)
        names = [n for n in extra.split(",") if n] + [n for n in FEATURES if n not in extra.split(",")]
        search.setMetaValue("extra_features", ",".join(names))
        prot.setSearchParameters(search)

    try:
        from pyopenms import PeptideIdentificationList
        container = PeptideIdentificationList()
        for pep in peptides:
            container.push_back(pep)
    except ImportError:
        container = peptides
    IdXMLFile().store(str(idxml_out), proteins, container)

    predicted = int(psms["predictable"].sum()) if not psms.empty else 0
    summary = {
        "file": Path(mzml).name, "instrument": instrument, "nce": nce, "psms": len(psms),
        "predicted": predicted,
        "median_angle_targets": float(values.loc[~psms["decoy"], "peptdeep:spectral_angle"].median()) if len(values) else None,
        "median_angle_decoys": float(values.loc[psms["decoy"], "peptdeep:spectral_angle"].median()) if len(values) and psms["decoy"].any() else None,
    }
    log(
        f"{summary['file']}: {instrument}, NCE {nce:g}; {predicted} of {len(psms)} PSMs predicted; "
        f"median spectral angle targets {summary['median_angle_targets'] or 0:.2f}, "
        f"decoys {summary['median_angle_decoys'] or 0:.2f}"
    )
    return summary
