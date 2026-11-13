"""Predicted fragment spectra (AlphaPeptDeep) to compare with a measured PSM spectrum.

Streamlit-free: the page caches the model and renders the result.
"""
from __future__ import annotations

import importlib.util
import os
import sys
import time
import numpy as np
import pandas as pd
import polars as pl
from pyopenms import AASequence, MSExperiment, MzMLFile, PeakFileOptions

# peptdeep's default model setting; used when the instrument is not recognised
DEFAULT_INSTRUMENT = "Lumos"
DEFAULT_NCE = 30.0
# Tries of the first-use model download, and the wait before the next one
_DOWNLOAD_ATTEMPTS = 4
_DOWNLOAD_RETRY_SECONDS = 15
FRAG_TYPES = ["b_z1", "y_z1", "b_z2", "y_z2"]
# Ion labels as the SequenceView writes them: y7¹⁺, b3²⁺
_SUPERSCRIPT = {1: "¹", 2: "²"}

# Substrings of mzML instrument names -> peptdeep instrument groups
_INSTRUMENT_HINTS = [
    ("astral", "Astral"),
    ("exploris", "Exploris"),
    ("q exactive", "QE"),
    ("qexactive", "QE"),
    ("lumos", "Lumos"),
    ("fusion", "Lumos"),
    ("eclipse", "Lumos"),
    ("ascend", "Lumos"),
    ("timstof", "timsTOF"),
    ("tims tof", "timsTOF"),
    ("tripletof", "SciexTOF"),
    ("zenotof", "SciexTOF"),
    ("velos", "Velos"),
    ("elite", "Elite"),
]


def peptdeep_available() -> bool:
    """Whether peptdeep is installed, without importing it.

    Importing peptdeep loads torch, which must stay out of the app server: on
    Windows, polars' parquet scan crashed the server (access violation) once
    torch was loaded in the same process. Predictions run in their own process
    (``prediction_worker``).
    """
    if "peptdeep" in sys.modules:
        return True
    try:
        return importlib.util.find_spec("peptdeep") is not None
    except (ImportError, ValueError):
        return False


def load_model():
    """The pretrained AlphaPeptDeep models on the CPU (downloaded once on first use)."""
    from peptdeep.pretrained_models import ModelManager, get_model_zip_file_path
    from peptdeep.settings import global_settings

    # The Docker images ship the models outside the (per-user) home directory
    if os.environ.get("PEPTDEEP_HOME"):
        global_settings["PEPTDEEP_HOME"] = os.environ["PEPTDEEP_HOME"]
    # The first use downloads the models from GitHub; a slow or failing
    # download raises FileNotFoundError, so try it a few times. A download cut
    # off while reading leaves a partial zip, which peptdeep would not fetch
    # again, so it goes first.
    for attempt in range(_DOWNLOAD_ATTEMPTS):
        try:
            model = ModelManager(device="cpu")
            break
        except FileNotFoundError:
            if attempt == _DOWNLOAD_ATTEMPTS - 1:
                raise
            try:
                os.remove(get_model_zip_file_path())
            except FileNotFoundError:
                pass
            time.sleep(_DOWNLOAD_RETRY_SECONDS * (attempt + 1))
    model.load_installed_models()
    return model


def instrument_of(mzml_path) -> str:
    """peptdeep instrument group for the instrument named in an mzML header."""
    options = PeakFileOptions()
    options.setMetadataOnly(True)
    mzml = MzMLFile()
    mzml.setOptions(options)
    exp = MSExperiment()
    try:
        mzml.load(str(mzml_path), exp)
    except Exception:
        return DEFAULT_INSTRUMENT
    name = exp.getInstrument().getName().lower()
    for hint, group in _INSTRUMENT_HINTS:
        if hint in name:
            return group
    return DEFAULT_INSTRUMENT


def _mod_name(name: str, site: str, n_term: bool, mod_index) -> str | None:
    candidates = [f"{name}@{site}"]
    if n_term:
        candidates.append(f"{name}@{site}^Any_N-term")
    for candidate in candidates:
        if candidate in mod_index:
            return candidate
    return None


def to_peptdeep(peptide: str) -> tuple[str, str, str] | None:
    """OpenMS peptide string -> (sequence, mods, mod_sites) in AlphaBase notation.

    None when a modification has no AlphaBase name, so it cannot be predicted.
    """
    from alphabase.constants.modification import MOD_DF

    seq = AASequence.fromString(peptide)
    mods, sites = [], []
    if seq.hasNTerminalModification():
        term, first = seq.getNTerminalModificationName(), seq.getResidue(0).getOneLetterCode()
        # residue-specific terminal mods (pyro-Glu) are named after their residue
        candidates = (f"{term}@Any_N-term", f"{term}@Protein_N-term",
                      f"{term}@{first}^Any_N-term", f"{term}@{first}^Protein_N-term")
        name = next((m for m in candidates if m in MOD_DF.index), None)
        if name is None:
            return None
        mods.append(name)
        sites.append("0")
    for i in range(seq.size()):
        residue = seq.getResidue(i)
        if not residue.isModified():
            continue
        name = _mod_name(residue.getModificationName(), residue.getOneLetterCode(), i == 0, MOD_DF.index)
        if name is None:
            return None
        mods.append(name)
        sites.append(str(i + 1))
    if seq.hasCTerminalModification():
        term = seq.getCTerminalModificationName()
        name = next((m for m in (f"{term}@Any_C-term", f"{term}@Protein_C-term") if m in MOD_DF.index), None)
        if name is None:
            return None
        mods.append(name)
        sites.append("-1")
    return seq.toUnmodifiedString(), ";".join(mods), ";".join(sites)


def predict_fragments(model, peptide: str, charge: int, instrument: str = DEFAULT_INSTRUMENT,
                      nce: float = DEFAULT_NCE) -> pd.DataFrame | None:
    """Predicted b/y ions of one precursor: columns mz, intensity (max 1), ion.

    Fragment charges stop below the precursor charge. None when the peptide
    cannot be predicted (unknown modification, unsupported length or residue).
    """
    converted = to_peptdeep(peptide)
    if converted is None:
        return None
    sequence, mods, sites = converted
    if not 3 <= len(sequence) <= 60 or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        return None
    precursor = pd.DataFrame({
        "sequence": [sequence], "mods": [mods], "mod_sites": [sites],
        "charge": [int(charge)], "nce": [float(nce)], "instrument": [instrument],
    })
    result = model.predict_all(precursor, predict_items=["ms2"], frag_types=FRAG_TYPES, multiprocessing=False)
    mzs, intensities = result["fragment_mz_df"], result["fragment_intensity_df"]
    n_sites = len(mzs)
    rows = []
    for frag_type in FRAG_TYPES:
        ion, z = frag_type.split("_z")
        z = int(z)
        if z > 1 and z >= int(charge):
            continue
        for i in range(n_sites):
            intensity = float(intensities[frag_type].iat[i])
            if intensity <= 0:
                continue
            # peptdeep rows run along the backbone: b ions grow, y ions shrink
            number = i + 1 if ion == "b" else n_sites - i
            rows.append({"mz": float(mzs[frag_type].iat[i]), "intensity": intensity,
                         "ion": f"{ion}{number}{_SUPERSCRIPT[z]}⁺"})
    if not rows:
        return None
    predicted = pd.DataFrame(rows).sort_values("mz", ignore_index=True)
    predicted["intensity"] /= predicted["intensity"].max()
    return predicted


def match_peaks(measured_mz: np.ndarray, measured_int: np.ndarray, predicted_mz: np.ndarray,
                tolerance: float, ppm: bool) -> np.ndarray:
    """For each predicted ion, the most intense measured intensity within tolerance (0 if none)."""
    order = np.argsort(measured_mz)
    mz, inten = np.asarray(measured_mz)[order], np.asarray(measured_int)[order]
    matched = np.zeros(len(predicted_mz))
    for k, target in enumerate(predicted_mz):
        tol = target * tolerance * 1e-6 if ppm else tolerance
        lo, hi = np.searchsorted(mz, target - tol), np.searchsorted(mz, target + tol, side="right")
        if hi > lo:
            matched[k] = inten[lo:hi].max()
    return matched


def spectral_angle(observed: np.ndarray, predicted: np.ndarray) -> float:
    """Normalised spectral contrast angle (Prosit): 1 identical, 0 nothing in common."""
    a, b = np.asarray(observed, float), np.asarray(predicted, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    cosine = float(np.clip(a @ b / (na * nb), -1.0, 1.0))
    return 1 - 2 * np.arccos(cosine) / np.pi


def cosine_similarity(observed: np.ndarray, predicted: np.ndarray) -> float:
    """Cosine of the measured and predicted fragment intensities: 1 identical, 0 nothing in common."""
    a, b = np.asarray(observed, float), np.asarray(predicted, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return float(np.clip(a @ b / (na * nb), 0.0, 1.0))


def mirror_frame(measured: pl.DataFrame, predicted: pd.DataFrame) -> pl.DataFrame:
    """One frame for MirrorPlot: measured peaks (side 0, their own peak_id) and
    predicted ions (side 1, negative peak_id), both scaled to 100 % of their base peak.
    """
    top = measured.select(
        pl.lit(0, pl.Int8).alias("side"),
        pl.col("peak_id").cast(pl.Int64),
        pl.col("mass").cast(pl.Float64),
        (pl.col("intensity") / pl.col("intensity").max() * 100).cast(pl.Float64).alias("intensity"),
        pl.lit(False).alias("highlight"),
        pl.lit("").alias("ion"),
    )
    bottom = pl.DataFrame({
        "side": pl.Series([1] * len(predicted), dtype=pl.Int8),
        "peak_id": pl.Series(-1 - np.arange(len(predicted)), dtype=pl.Int64),
        "mass": predicted["mz"].to_numpy(float),
        "intensity": predicted["intensity"].to_numpy(float) * 100,
        "highlight": [True] * len(predicted),
        "ion": predicted["ion"].astype(str).tolist(),
    })
    return pl.concat([top, bottom])
