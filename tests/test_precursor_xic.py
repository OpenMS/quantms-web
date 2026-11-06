"""Precursor XICs of a PSM (src/common/precursor_xic.py), on a small synthetic mzML."""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

oms = pytest.importorskip("pyopenms")

from src.common import precursor_xic as xic  # noqa: E402

MZ, CHARGE = 500.0, 2


def _write_mzml(path, apex=100.0, scale=1.0, extra=None):
    """MS1 spectra every 2 s from 0 to 200 s: a Gaussian elution of MZ and its isotopes, plus an MS2."""
    exp = oms.MSExperiment()
    for rt in np.arange(0.0, 201.0, 2.0):
        s = oms.MSSpectrum()
        s.setMSLevel(1)
        s.setRT(float(rt))
        height = 1e6 * scale * np.exp(-((rt - apex) ** 2) / (2 * 8.0**2))
        mzs = [MZ + k * xic.C13_DELTA / CHARGE for k in range(3)] + [MZ + 0.2, 650.0]
        ints = [height, 0.6 * height, 0.2 * height, 5e5, 7e5]
        if extra is not None and rt == extra:
            ints[0] += 3e6
        s.set_peaks((np.array(mzs), np.array(ints)))
        exp.addSpectrum(s)
        if rt == apex:
            ms2 = oms.MSSpectrum()
            ms2.setMSLevel(2)
            ms2.setRT(float(rt) + 0.5)
            ms2.set_peaks((np.array([MZ]), np.array([1e9])))
            exp.addSpectrum(ms2)
    oms.MzMLFile().store(str(path), exp)
    return path


def test_isotope_traces(tmp_path):
    mzml = _write_mzml(tmp_path / "a.mzML")
    traces = xic.extract_xic(mzml, MZ, CHARGE, rt=101.0)
    assert set(traces["series"]) == {"M", "M+1", "M+2"}
    mono = traces[traces["series"] == "M"]
    # ±60 s around 101 s, one point per MS1 spectrum; the MS2 is not read
    assert mono["rt"].min() >= 41.0 and mono["rt"].max() <= 161.0
    assert len(mono) == 60
    # the neighbouring peak 0.2 m/z away stays outside 10 ppm
    assert mono["intensity"].max() == pytest.approx(1e6)
    summary = xic.summarize(traces, ms2_rt=100.5)
    assert summary["apex_rt"] == 100.0
    assert summary["ms2_offset"] == pytest.approx(0.5)
    assert summary["isotope_correlation"] == pytest.approx(1.0)
    assert summary["area"] > 0


def test_interference_lowers_correlation(tmp_path):
    mzml = _write_mzml(tmp_path / "b.mzML", extra=120.0)
    summary = xic.summarize(xic.extract_xic(mzml, MZ, CHARGE, rt=100.0), ms2_rt=100.0)
    assert summary["apex_rt"] == 120.0
    assert summary["isotope_correlation"] < 0.8  # a clean peak gives 1.0


def test_across_runs_and_marker(tmp_path):
    runs = {
        "ctrl": _write_mzml(tmp_path / "ctrl.mzML"),
        "drug": _write_mzml(tmp_path / "drug.mzML", apex=130.0, scale=4.0),
        "missing": tmp_path / "missing.mzML",
    }
    traces = xic.xic_across_runs(runs, MZ, CHARGE, rt=100.0)
    assert set(traces["series"]) == {"ctrl", "drug"}
    peaks = traces.groupby("series")["intensity"].max()
    assert peaks["drug"] == pytest.approx(4 * peaks["ctrl"])
    marked = xic.with_ms2_marker(traces, 100.0)
    line = marked[marked["series"] == xic.MS2_SERIES]
    assert line["rt"].tolist() == [100.0, 100.0]
    assert line["intensity"].max() == pytest.approx(peaks["drug"])


def test_precursor_mz_from_sequence():
    mz = xic.precursor_mz("PEPTIDEK", 2)
    assert mz == pytest.approx(oms.AASequence.fromString("PEPTIDEK").getMZ(2))
    assert xic.precursor_mz("not a peptide ][", 2, fallback=123.4) == 123.4


def test_empty_trace():
    summary = xic.summarize(pd.DataFrame(columns=xic.XIC_COLUMNS), 10.0)
    assert summary["apex_rt"] is None and summary["area"] == 0.0
