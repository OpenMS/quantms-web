"""The MS1 peak map: precursors of PSMs marked on their MS1 peaks, and the chromatograms."""
from pathlib import Path

import polars as pl
import pytest

from src.common import peak_map


def _peaks():
    # Two MS1 spectra at 10 s and 20 s
    return pl.DataFrame({
        "rt": [10.0, 10.0, 10.0, 20.0, 20.0],
        "mz": [500.0, 500.002, 600.0, 500.001, 700.0],
        "intensity": [100.0, 300.0, 50.0, 200.0, 10.0],
        "ms1": [0, 0, 0, 1, 1],
    }, schema=peak_map.PEAK_SCHEMA)


def test_precursor_is_the_most_intense_peak_in_the_ms1_spectrum_before_the_psm():
    psms = pl.DataFrame({"id_idx": [0, 1], "mz": [500.001, 600.0], "rt": [15.0, 12.0]})
    found = peak_map.precursor_peaks(_peaks(), psms)
    # PSM 0 (MS2 at 15 s) was picked from the spectrum at 10 s: of its two peaks within 10 ppm the higher one
    assert dict(zip(found["id_idx"], found["row"])) == {0: 1, 1: 2}


def test_psm_after_the_second_spectrum_uses_it_and_no_peak_within_tolerance_is_left_out():
    psms = pl.DataFrame({"id_idx": [3, 4, 5], "mz": [500.0, 650.0, 500.0], "rt": [25.0, 25.0, 5.0]})
    found = peak_map.precursor_peaks(_peaks(), psms)
    # 650 has no peak; the PSM at 5 s comes before any MS1 spectrum
    assert found.to_dicts() == [{"id_idx": 3, "row": 3}]


def test_mark_precursors_keeps_every_peak_and_the_first_psm_of_a_shared_peak():
    psms = pl.DataFrame({"id_idx": [7, 8], "mz": [700.0, 700.0], "rt": [21.0, 21.0]})
    marked = peak_map.mark_precursors(_peaks(), psms)
    assert marked.height == 5
    assert marked["id_idx"].to_list() == [None, None, None, None, 7]
    assert marked.schema["id_idx"] == pl.Int64


def test_no_psms_marks_nothing():
    marked = peak_map.mark_precursors(_peaks(), pl.DataFrame(schema={"id_idx": pl.Int64, "mz": pl.Float64,
                                                                     "rt": pl.Float64}))
    assert marked["id_idx"].null_count() == 5


def test_chromatograms_sum_and_maximum_per_spectrum():
    chrom = peak_map.chromatograms(_peaks())
    tic = chrom.filter(pl.col("series") == peak_map.TIC).sort("rt")
    bpc = chrom.filter(pl.col("series") == peak_map.BPC).sort("rt")
    assert tic["intensity"].to_list() == [450.0, 210.0]
    assert bpc["intensity"].to_list() == [300.0, 200.0]
    assert tic["rt"].to_list() == [10.0, 20.0]


def test_ms1_peaks_reads_only_ms1(tmp_path: Path):
    oms = pytest.importorskip("pyopenms")
    exp = oms.MSExperiment()
    for rt, level, mzs in [(1.0, 1, [400.0, 401.0]), (2.0, 2, [150.0]), (3.0, 1, [402.0])]:
        spectrum = oms.MSSpectrum()
        spectrum.setRT(rt)
        spectrum.setMSLevel(level)
        spectrum.set_peaks((mzs, [10.0] * len(mzs)))
        exp.addSpectrum(spectrum)
    path = tmp_path / "run.mzML"
    oms.MzMLFile().store(str(path), exp)
    peaks = peak_map.ms1_peaks(path)
    assert peaks["mz"].to_list() == pytest.approx([400.0, 401.0, 402.0])
    assert peaks["ms1"].to_list() == [0, 0, 1]
