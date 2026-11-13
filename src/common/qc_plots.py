"""Identification QC: score distribution, precursor mass error and PSMs per file.

``psm_qc_table`` reads one idXML into a flat pandas table (one row per PSM, top
hit only); the plot functions take the concatenated tables of all runs and
return Plotly figures, so the page only decides where to show them.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from src.common.results_helpers import extract_filename_from_idxml, load_idxml

# Spacing of the 13C isotope peaks (Da): the precursor may have been picked on one of them
C13_C12_DIFF = 1.0033548
MAX_ISOTOPE_ERROR = 3

# Percolator writes its scores as PSI-MS meta values on each hit
PERCOLATOR_SCORE_KEYS = ("MS:1001492", "Percolator_score")
Q_VALUE_KEYS = ("MS:1001491", "q-value", "Percolator_qvalue")

TARGET_COLOR = "#2a78d6"
DECOY_COLOR = "#eb6834"

QC_COLUMNS = ["file", "decoy", "score", "score_type", "percolator_score", "q_value", "charge", "precursor_error_ppm"]


def _meta(obj, keys):
    for key in keys:
        if obj.metaValueExists(key):
            try:
                return float(obj.getMetaValue(key))
            except (TypeError, ValueError):
                return np.nan
    return np.nan


def _is_decoy(hit, decoy_prefixes=("DECOY_", "rev_")) -> bool:
    """A hit is a decoy if PeptideIndexer said so, else if every protein carries a decoy prefix."""
    if hit.metaValueExists("target_decoy"):
        value = hit.getMetaValue("target_decoy")
        if isinstance(value, bytes):
            value = value.decode()
        return str(value) == "decoy"
    accessions = [ev.getProteinAccession() for ev in hit.getPeptideEvidences()]
    return bool(accessions) and all(a.startswith(decoy_prefixes) for a in accessions)


def precursor_error_ppm(observed_mz: float, theoretical_mz: float, charge: int) -> float:
    """Precursor error in ppm, after removing a 13C isotope error of up to three peaks."""
    if not charge or theoretical_mz <= 0:
        return np.nan
    isotope = round((observed_mz - theoretical_mz) * charge / C13_C12_DIFF)
    if abs(isotope) <= MAX_ISOTOPE_ERROR:
        observed_mz -= isotope * C13_C12_DIFF / charge
    return (observed_mz - theoretical_mz) / theoretical_mz * 1e6


def psm_qc_table(idxml_path: Path, mzml_names=None) -> pd.DataFrame:
    """One row per spectrum (its top hit) of an idXML file."""
    _, peptides = load_idxml(idxml_path)
    file = extract_filename_from_idxml(Path(idxml_path), mzml_names)
    records = []
    for pep in peptides:
        hits = pep.getHits()
        if not hits:
            continue
        hit = hits[0]
        charge = hit.getCharge()
        theoretical = hit.getSequence().getMZ(charge) if charge else np.nan
        records.append({
            "file": file,
            "decoy": _is_decoy(hit),
            "score": hit.getScore(),
            "score_type": pep.getScoreType(),
            "percolator_score": _meta(hit, PERCOLATOR_SCORE_KEYS),
            "q_value": _meta(hit, Q_VALUE_KEYS),
            "charge": charge,
            "precursor_error_ppm": precursor_error_ppm(pep.getMZ(), theoretical, charge),
        })
    return pd.DataFrame(records, columns=QC_COLUMNS)


def _score_axis(df: pd.DataFrame) -> tuple[pd.Series, str]:
    """The score to plot, higher is better: Percolator's SVM score, else the search engine's."""
    if df["percolator_score"].notna().all():
        return df["percolator_score"], "Percolator score"
    score_type = str(df["score_type"].iloc[0]) if len(df) else "score"
    if "e-value" in score_type.lower() or "expect" in score_type.lower():
        return -np.log10(df["score"].clip(lower=1e-300)), "-log10(E-value)"
    return df["score"], score_type


def score_distribution(df: pd.DataFrame, fdr: float | None = None) -> go.Figure:
    """Overlaid histograms of target and decoy PSM scores.

    Decoys should sit at low scores and targets show a second mode of correct
    matches at high scores; the dashed line marks the lowest target score that
    passes the PSM FDR.
    """
    scores, label = _score_axis(df)
    fig = go.Figure()
    for is_decoy, name, color in ((False, "Targets", TARGET_COLOR), (True, "Decoys", DECOY_COLOR)):
        values = scores[df["decoy"] == is_decoy].dropna()
        fig.add_trace(go.Histogram(x=values, name=f"{name} ({len(values):,})", marker_color=color, opacity=0.6, nbinsx=60))
    if fdr is not None and fdr < 1.0:
        passing = scores[(~df["decoy"]) & (df["q_value"] <= fdr)].dropna()
        if len(passing):
            fig.add_vline(
                x=passing.min(), line_dash="dash", line_color="grey",
                annotation_text=f"PSM FDR {fdr * 100:.3g} % (q ≤ {fdr:.3g})",
            )
    fig.update_layout(
        barmode="overlay",
        title="Score distribution of target and decoy PSMs",
        xaxis_title=label,
        yaxis_title="PSMs",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
    )
    return fig


def precursor_error_histogram(df: pd.DataFrame) -> go.Figure:
    """Precursor mass error of target PSMs; the median shows the calibration offset."""
    errors = df.loc[~df["decoy"], "precursor_error_ppm"].dropna()
    fig = go.Figure(go.Histogram(x=errors, marker_color=TARGET_COLOR, nbinsx=60, name="Target PSMs"))
    if len(errors):
        median = float(errors.median())
        fig.add_vline(x=median, line_dash="dash", line_color="grey", annotation_text=f"median {median:+.2f} ppm")
        radius = max(1.0, float(errors.abs().quantile(0.99)))
        fig.update_xaxes(range=[-radius * 1.05, radius * 1.05])
    fig.update_layout(
        title=f"Precursor mass error of target PSMs ({len(errors):,})",
        xaxis_title="Precursor error (ppm), isotope errors removed",
        yaxis_title="PSMs",
        showlegend=False,
    )
    return fig


def psms_per_file(df: pd.DataFrame) -> pd.DataFrame:
    """Target PSMs per run (file, PSMs), so a failed or weak run stands out."""
    files = sorted(df["file"].unique())
    counts = df.loc[~df["decoy"]].groupby("file").size().reindex(files, fill_value=0)
    return pd.DataFrame({"file": files, "PSMs": counts.to_numpy()})


def _binned(values: pd.Series, edges: np.ndarray) -> np.ndarray:
    counts, _ = np.histogram(values.dropna().to_numpy(dtype=float), bins=edges)
    return counts


def score_histogram(df: pd.DataFrame, bins: int = 60) -> tuple[pd.DataFrame, str, float | None]:
    """Target and decoy PSMs per score bin as a tidy table (x, y, kind) for a line plot.

    Returns the table, the score's name and None (the FDR cut-off is added by
    ``fdr_cutoff``).
    """
    scores, label = _score_axis(df)
    finite = scores.replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return pd.DataFrame(columns=["x", "y", "kind"]), label, None
    edges = np.linspace(float(finite.min()), float(finite.max()) or 1.0, bins + 1)
    if edges[0] == edges[-1]:
        edges = np.linspace(edges[0] - 0.5, edges[0] + 0.5, bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    frames = [
        pd.DataFrame({"x": centers, "y": _binned(scores[df["decoy"] == decoy], edges).astype(float), "kind": kind})
        for decoy, kind in ((False, "target"), (True, "decoy"))
    ]
    return pd.concat(frames, ignore_index=True), label, None


def fdr_cutoff(df: pd.DataFrame, fdr: float) -> float | None:
    """The lowest target score that passes the PSM FDR (q-value ≤ fdr)."""
    scores, _ = _score_axis(df)
    passing = scores[(~df["decoy"]) & (df["q_value"] <= fdr)].dropna()
    return float(passing.min()) if len(passing) and fdr < 1.0 else None


def mass_error_histogram(df: pd.DataFrame, bins: int = 60) -> tuple[pd.DataFrame, float | None]:
    """Target PSMs per precursor error bin (x, y) and the median error, for a line plot."""
    errors = df.loc[~df["decoy"], "precursor_error_ppm"].dropna()
    if errors.empty:
        return pd.DataFrame(columns=["x", "y"]), None
    radius = max(1.0, float(errors.abs().quantile(0.99))) * 1.05
    edges = np.linspace(-radius, radius, bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    return pd.DataFrame({"x": centers, "y": _binned(errors.clip(-radius, radius), edges).astype(float)}), float(errors.median())
