"""build_spectra_cache and peaks_for_identifications."""
from pathlib import Path

import pyopenms as oms

import polars as pl

from src.common.results_helpers import build_spectra_cache, peaks_for_identifications


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


def test_peaks_for_identifications_keeps_identified_spectra_only():
    spectra = pl.DataFrame({
        "peak_id": [0, 1, 2, 3, 4],
        "file_index": [0, 0, 0, 1, 1],
        "scan_id": [2, 2, 3, 2, 3],
        "mass": [100.0, 200.0, 150.0, 100.0, 150.0],
        "intensity": [1.0, 2.0, 5.0, 1.0, 5.0],
    })
    ids = pl.DataFrame({"file_index": [0, 0, 1], "scan_id": [2, 2, 3]})

    out = peaks_for_identifications(spectra, ids).collect()

    assert out["peak_id"].to_list() == [0, 1, 4]


def test_write_spectra_cache_one_parquet_per_file(tmp_path):
    from src.common.results_helpers import write_spectra_cache

    a, b = tmp_path / "a.mzML", tmp_path / "b.mzML"
    _write_mzml(a)
    _write_mzml(b)
    index = {"a.mzML": 0, "b.mzML": 1}

    cache = write_spectra_cache(tmp_path / "cache", index, [a, b])

    assert sorted(cache) == [0, 1]
    whole, _ = build_spectra_cache(tmp_path, dict(index), [a, b])
    per_file = pl.concat([pl.read_parquet(cache[i]) for i in (0, 1)])
    assert per_file.equals(whole)  # same peaks and peak_ids as the in-memory cache
    assert write_spectra_cache(tmp_path / "empty", {}, []) == {}


def test_peaks_for_identifications_reads_only_the_identified_files(tmp_path):
    from src.common.results_helpers import write_spectra_cache

    a = tmp_path / "a.mzML"
    _write_mzml(a)
    cache = write_spectra_cache(tmp_path / "cache", {"a.mzML": 0}, [a])
    cache[1] = tmp_path / "missing.parquet"  # never read: no identification in file 1

    ids = pl.DataFrame({"file_index": [0], "scan_id": [3]})
    out = peaks_for_identifications(cache, ids).collect()

    assert out["mass"].to_list() == [150.0]
    none = peaks_for_identifications(cache, pl.DataFrame({"file_index": [5], "scan_id": [3]}, schema={"file_index": pl.Int64, "scan_id": pl.Int64}))
    assert none.collect().height == 0
