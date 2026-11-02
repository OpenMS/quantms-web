"""build_spectra_cache: one row per MS2 peak, MS1 skipped, columns typed."""
from pathlib import Path

import pyopenms as oms

from src.common.results_helpers import build_spectra_cache


def _write_mzml(path: Path) -> None:
    exp = oms.MSExperiment()
    ms1 = oms.MSSpectrum()
    ms1.setMSLevel(1)
    ms1.setNativeID("scan=1")
    ms1.set_peaks(([400.0, 500.0, 600.0], [10.0, 20.0, 30.0]))
    exp.addSpectrum(ms1)
    for scan, peaks in ((2, ([100.0, 200.0], [1.0, 2.0])), (3, ([150.0], [5.0])), (4, ([], []))):
        ms2 = oms.MSSpectrum()
        ms2.setMSLevel(2)
        ms2.setNativeID(f"scan={scan}")
        ms2.set_peaks(peaks)
        exp.addSpectrum(ms2)
    oms.MzMLFile().store(str(path), exp)


def test_ms2_peaks_only(tmp_path):
    a, b = tmp_path / "a.mzML", tmp_path / "b.mzML"
    _write_mzml(a)
    _write_mzml(b)

    df, index = build_spectra_cache(tmp_path, {"a.mzML": 0}, [a, b])

    assert index == {"a.mzML": 0, "b.mzML": 1}
    assert df.columns == ["peak_id", "file_index", "scan_id", "mass", "intensity"]
    assert df.height == 6  # 3 MS2 peaks per file, MS1 and the empty MS2 skipped
    assert df["peak_id"].to_list() == list(range(6))
    assert df["file_index"].to_list() == [0, 0, 0, 1, 1, 1]
    assert df["scan_id"].to_list() == [2, 2, 3, 2, 2, 3]
    assert df["mass"].to_list() == [100.0, 200.0, 150.0] * 2
    assert df["intensity"].to_list() == [1.0, 2.0, 5.0] * 2


def test_no_files(tmp_path):
    df, index = build_spectra_cache(tmp_path, {}, [])
    assert df.height == 0
    assert index == {}
