"""Extracted ion chromatograms (XICs) of a PSM's precursor.

Streamlit-free. The precursor of a peptide-spectrum match is traced over the
MS1 spectra of its run: the monoisotopic peak and the next isotope peaks, each
summed within a ppm window, over a retention time window around the MS2
spectrum. A good quantification shows one smooth elution peak whose isotope
traces rise and fall together; a second peak, a shoulder or isotope traces
that do not follow the monoisotopic one point to an interference.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

C13_DELTA = 1.0033548378
PROTON = 1.007276466812
ISOTOPES = 3
PPM = 10.0
RT_WINDOW = 60.0  # seconds on either side of the MS2 spectrum
RUNS_RT_WINDOW = 180.0  # wider across runs: their retention times are not aligned
MS2_SERIES = "MS2 spectrum"
_trapezoid = getattr(np, "trapezoid", None) or np.trapz
XIC_COLUMNS = ["rt", "intensity", "series"]


def isotope_label(k: int) -> str:
    """Series name of the k-th isotope trace: ``M``, ``M+1``, ..."""
    return "M" if k == 0 else f"M+{k}"


def precursor_mz(sequence: str, charge: int, fallback: float | None = None) -> float | None:
    """Monoisotopic m/z of a peptide at ``charge``.

    Taken from the sequence, so a precursor the instrument picked on an isotope
    peak is still traced from its monoisotopic peak. ``fallback`` (the m/z the
    search recorded) when the sequence cannot be read.
    """
    try:
        import pyopenms as oms  # pylint: disable=import-outside-toplevel

        return float(oms.AASequence.fromString(str(sequence)).getMZ(int(charge)))
    except Exception:  # pylint: disable=broad-except
        return fallback


def _ms1_spectra(mzml: Path, rt_low: float, rt_high: float, mz_low: float, mz_high: float):
    """The MS1 spectra of ``mzml`` between two retention times, cut to an m/z range."""
    import pyopenms as oms  # pylint: disable=import-outside-toplevel

    def span(low: float, high: float):
        # pyOpenMS 3.6 (nanobind) takes the bounds, earlier versions DPosition1s
        try:
            return oms.DRange1(float(low), float(high))
        except TypeError:
            return oms.DRange1(oms.DPosition1(low), oms.DPosition1(high))

    exp = oms.MSExperiment()
    reader = oms.MzMLFile()
    options = reader.getOptions()
    options.setMSLevels([1])
    options.setRTRange(span(rt_low, rt_high))
    options.setMZRange(span(mz_low, mz_high))
    options.setFillData(True)
    reader.setOptions(options)
    reader.load(str(mzml), exp)
    return [s for s in exp.getSpectra() if s.getMSLevel() == 1]


def extract_xic(
    mzml: Path | str,
    mz: float,
    charge: int,
    rt: float,
    *,
    ppm: float = PPM,
    rt_window: float = RT_WINDOW,
    isotopes: int = ISOTOPES,
) -> pd.DataFrame:
    """Intensity of the precursor's isotope peaks in every MS1 spectrum near ``rt``.

    One row per MS1 spectrum and isotope (``rt``, ``intensity``, ``series``
    as named by ``isotope_label``); spectra without the peak count as 0, so
    a trace that drops out shows as a gap to the baseline.
    """
    charge = max(int(charge), 1)
    targets = [mz + k * C13_DELTA / charge for k in range(isotopes)]
    margin = targets[-1] * ppm * 1e-6
    spectra = _ms1_spectra(Path(mzml), rt - rt_window, rt + rt_window, targets[0] - margin, targets[-1] + margin)
    rows = []
    for spectrum in spectra:
        mzs, intensities = spectrum.get_peaks()
        for k, target in enumerate(targets):
            tolerance = target * ppm * 1e-6
            inside = np.abs(np.asarray(mzs) - target) <= tolerance
            rows.append((float(spectrum.getRT()), float(np.asarray(intensities)[inside].sum()), isotope_label(k)))
    return pd.DataFrame(rows, columns=XIC_COLUMNS).sort_values(["series", "rt"], kind="stable", ignore_index=True)


def xic_across_runs(
    mzmls: dict[str, Path | str],
    mz: float,
    charge: int,
    rt: float,
    *,
    ppm: float = PPM,
    rt_window: float = RUNS_RT_WINDOW,
) -> pd.DataFrame:
    """The monoisotopic trace of the precursor in every run, named by run (``series``)."""
    frames = []
    for name, path in mzmls.items():
        if path is None or not Path(path).is_file():
            continue
        trace = extract_xic(path, mz, charge, rt, ppm=ppm, rt_window=rt_window, isotopes=1)
        frames.append(trace.assign(series=name))
    if not frames:
        return pd.DataFrame(columns=XIC_COLUMNS)
    return pd.concat(frames, ignore_index=True)


def with_ms2_marker(xic: pd.DataFrame, ms2_rt: float) -> pd.DataFrame:
    """``xic`` plus a vertical line at the retention time of the MS2 spectrum."""
    top = float(xic["intensity"].max()) if len(xic) else 0.0
    marker = pd.DataFrame(
        {"rt": [ms2_rt, ms2_rt], "intensity": [0.0, top or 1.0], "series": [MS2_SERIES, MS2_SERIES]}
    )
    return pd.concat([xic, marker], ignore_index=True)


def summarize(xic: pd.DataFrame, ms2_rt: float) -> dict:
    """Apex and area of the monoisotopic trace, and where the MS2 spectrum fell on it.

    ``isotope_correlation`` is the Pearson correlation of the M+1 trace with
    the M trace (``None`` with too few points): close to 1 for a clean peak.
    """
    mono = xic[xic["series"] == isotope_label(0)].sort_values("rt")
    if mono.empty or not float(mono["intensity"].max()):
        return {"apex_rt": None, "apex_intensity": 0.0, "area": 0.0, "ms2_offset": None, "isotope_correlation": None}
    apex = mono.loc[mono["intensity"].idxmax()]
    area = float(_trapezoid(mono["intensity"].to_numpy(), mono["rt"].to_numpy())) if len(mono) > 1 else 0.0
    correlation = None
    second = xic[xic["series"] == isotope_label(1)].sort_values("rt")
    if len(second) == len(mono) and len(mono) >= 3:
        a, b = mono["intensity"].to_numpy(), second["intensity"].to_numpy()
        if a.std() > 0 and b.std() > 0:
            correlation = float(np.corrcoef(a, b)[0, 1])
    return {
        "apex_rt": float(apex["rt"]),
        "apex_intensity": float(apex["intensity"]),
        "area": area,
        "ms2_offset": float(ms2_rt) - float(apex["rt"]),
        "isotope_correlation": correlation,
    }
