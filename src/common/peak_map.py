"""The MS1 peak map of a run, with the precursors of its PSMs marked, and its chromatograms.

Streamlit-free. ``ms1_peaks`` reads every MS1 peak of an mzML file (retention
time, m/z, intensity, the MS1 spectrum it belongs to). ``mark_precursors``
gives the peak each PSM's MS2 spectrum was triggered on the PSM's row id
(``id_idx``): the most intense peak within a ppm window of the recorded
precursor m/z in the last MS1 spectrum before the MS2 spectrum. A click on
that peak in the map then picks the PSM, and a PSM picked elsewhere rings its
peak. ``chromatograms`` sums each MS1 spectrum (total ion current) and takes
its highest peak (base peak).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

PPM = 10.0
TIC, BPC = "Total ion current", "Base peak"
PEAK_SCHEMA = {"rt": pl.Float64, "mz": pl.Float64, "intensity": pl.Float64, "ms1": pl.Int64}
PRECURSOR_COLUMN = "id_idx"


def ms1_peaks(mzml: Path | str) -> pl.DataFrame:
    """Every peak of the MS1 spectra of ``mzml``: rt (s), mz, intensity, ms1 (spectrum number from 0)."""
    import pyopenms as oms  # pylint: disable=import-outside-toplevel

    exp = oms.MSExperiment()
    reader = oms.MzMLFile()
    options = reader.getOptions()
    options.setMSLevels([1])
    reader.setOptions(options)
    reader.load(str(mzml), exp)
    rts, mzs, intensities, numbers = [], [], [], []
    n = 0
    for spectrum in exp.getSpectra():
        if spectrum.getMSLevel() != 1:
            continue
        mz, intensity = spectrum.get_peaks()
        if len(mz):
            rts.append(np.full(len(mz), spectrum.getRT()))
            mzs.append(np.asarray(mz, dtype=float))
            intensities.append(np.asarray(intensity, dtype=float))
            numbers.append(np.full(len(mz), n, dtype=np.int64))
        n += 1
    if not rts:
        return pl.DataFrame(schema=PEAK_SCHEMA)
    peaks = pl.DataFrame({
        "rt": np.concatenate(rts), "mz": np.concatenate(mzs),
        "intensity": np.concatenate(intensities), "ms1": np.concatenate(numbers),
    }, schema=PEAK_SCHEMA)
    return peaks.filter(pl.col("intensity") > 0)


def precursor_peaks(peaks: pl.DataFrame, psms: pl.DataFrame, ppm: float = PPM) -> pl.DataFrame:
    """The MS1 peak each PSM was picked on: ``id_idx`` and the peak's row in ``peaks``.

    Args:
        peaks: ``ms1_peaks``.
        psms: ``id_idx``, ``mz`` (recorded precursor m/z) and ``rt`` (of the MS2 spectrum).

    PSMs without a peak within ``ppm`` in the MS1 spectrum before them are left out.
    """
    empty = pl.DataFrame(schema={PRECURSOR_COLUMN: pl.Int64, "row": pl.Int64})
    if peaks.is_empty() or psms.is_empty():
        return empty
    peaks = peaks.with_row_index("row")
    scans = peaks.group_by("ms1").agg(pl.col("rt").first()).sort("ms1")
    scan_rts = scans["rt"].to_numpy()
    scan_ids = scans["ms1"].to_numpy()
    by_scan = {int(k[0]): g for k, g in peaks.group_by("ms1")}
    found_ids, found_rows = [], []
    for psm_id, mz, rt in psms.select(PRECURSOR_COLUMN, "mz", "rt").iter_rows():
        if mz is None or rt is None:
            continue
        at = int(np.searchsorted(scan_rts, float(rt), side="right")) - 1
        if at < 0:
            continue
        scan = by_scan.get(int(scan_ids[at]))
        tolerance = float(mz) * ppm * 1e-6
        near = scan.filter((pl.col("mz") - float(mz)).abs() <= tolerance)
        if near.is_empty():
            continue
        found_ids.append(int(psm_id))
        found_rows.append(int(near.sort("intensity", descending=True)["row"][0]))
    if not found_ids:
        return empty
    return pl.DataFrame({PRECURSOR_COLUMN: found_ids, "row": found_rows}, schema=empty.schema)


def mark_precursors(peaks: pl.DataFrame, psms: pl.DataFrame, ppm: float = PPM) -> pl.DataFrame:
    """``peaks`` with ``id_idx`` set on the peak each PSM was picked on (null elsewhere).

    A peak several PSMs were picked on keeps the first of them (the best
    hit of the MS2 spectrum comes first in the idXML).
    """
    found = precursor_peaks(peaks, psms, ppm).unique("row", keep="first", maintain_order=True)
    marked = peaks.with_row_index("row").join(found, on="row", how="left").drop("row")
    return marked.with_columns(pl.col(PRECURSOR_COLUMN).cast(pl.Int64))


def chromatograms(peaks: pl.DataFrame) -> pl.DataFrame:
    """Total ion current and base peak per MS1 spectrum: rt, intensity, series."""
    if peaks.is_empty():
        return pl.DataFrame(schema={"rt": pl.Float64, "intensity": pl.Float64, "series": pl.Utf8})
    per_scan = peaks.group_by("ms1").agg(
        pl.col("rt").first(), pl.col("intensity").sum().alias(TIC), pl.col("intensity").max().alias(BPC),
    ).sort("ms1")
    return per_scan.unpivot(index="rt", on=[TIC, BPC], variable_name="series", value_name="intensity").select(
        "rt", "intensity", "series"
    )
