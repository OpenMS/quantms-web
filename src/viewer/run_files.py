"""One LC-MS run (an mzML file) of an experiment, for the Viewer's Single run.

An experiment is shown as a whole (Comparison) or one of its runs on its own.
The runs are the mzML files of every workflow run of the experiment, named in
the Viewer as ``<workflow run id>/<file stem>``. This module lists them, finds
each one's spectra and turns them into the columns the library's MS data
components read (``src/viz/ms_data``).
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

SEP = "/"


def file_id(run_id: str, stem: str) -> str:
    return f"{run_id}{SEP}{stem}"


def split(dataset: str) -> tuple[str, str]:
    """The workflow run and the file stem of a Single run dataset."""
    run_id, _, stem = dataset.partition(SEP)
    return run_id, stem


def is_file(dataset: str | None) -> bool:
    return bool(dataset) and SEP in dataset


def stems(params: dict) -> list[str]:
    """The file stems of a workflow run, in the order it processed them."""
    names = params.get("mzML-files") or []
    if isinstance(names, str):
        names = [names]
    return [Path(n).stem for n in names]


def mzml_path(params: dict, results_dir: Path, stem: str) -> tuple[Path | None, bool]:
    """Where the spectra of ``stem`` are read from, and whether it is the whole
    file (False: only the identified spectra the run kept)."""
    from src.common.identified_spectra import SPECTRA_DIR
    from src.common.spectrum_view import psm_mzml_files

    for path in psm_mzml_files(params, results_dir):
        if path is not None and path.stem == stem:
            return path, Path(path).parent != Path(results_dir, SPECTRA_DIR)
    return None, False


def read_tables(path: Path) -> dict[str, pl.DataFrame]:
    """The tables of one mzML file: ``scans`` (scan table), ``peaks`` (every
    peak by spectrum), ``ms1`` (the MS1 peak map) and ``chromatogram`` (TIC and
    base peak of the MS1 spectra)."""
    import numpy as np
    import pyopenms as oms

    from src.common.results_helpers import extract_scan_number

    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(path), exp)
    scans, peaks_mz, peaks_int, peaks_idx, chrom = [], [], [], [], []
    for index, spec in enumerate(exp.getSpectra()):
        mz, intensity = spec.get_peaks()
        level, rt = spec.getMSLevel(), spec.getRT()
        precursors = spec.getPrecursors()
        precursor = float("nan")
        if precursors:
            p = precursors[0]
            charge = p.getCharge()
            precursor = (p.getMZ() - 1.007276) * charge if charge else p.getMZ()
        scans.append({
            "index": index, "Scan": extract_scan_number(spec.getNativeID()) or index, "MSLevel": level,
            "RT": rt, "PrecursorMass": precursor, "#Peaks": len(mz),
        })
        peaks_mz.append(mz)
        peaks_int.append(intensity)
        peaks_idx.append(np.full(len(mz), index, dtype=np.int64))
        if level == 1:
            chrom.append((rt, float(intensity.sum()) if len(intensity) else 0.0,
                          float(intensity.max()) if len(intensity) else 0.0, index))
    levels = {s["index"]: s["MSLevel"] for s in scans}
    rts = {s["index"]: s["RT"] for s in scans}
    concat = (lambda parts, dtype: np.concatenate(parts).astype(dtype) if parts else np.array([], dtype=dtype))
    peaks = pl.DataFrame({
        "index": concat(peaks_idx, np.int64), "mz": concat(peaks_mz, np.float64),
        "intensity": concat(peaks_int, np.float64),
    })
    ms1 = (peaks.with_columns(pl.col("index").replace_strict(levels, default=0).alias("level"))
           .filter(pl.col("level") == 1)
           .select(pl.col("index").replace_strict(rts, default=0.0).alias("rt"),
                   pl.col("mz").alias("mass"), "intensity", pl.col("index").alias("scan_idx")))
    chromatogram = pl.DataFrame(
        [{"x": rt, "y": y, "series": name, "index": i}
         for rt, tic, bpc, i in chrom for name, y in (("TIC", tic), ("Base peak", bpc))],
        schema={"x": pl.Float64, "y": pl.Float64, "series": pl.Utf8, "index": pl.Int64},
    )
    return {"scans": pl.DataFrame(scans, schema={
        "index": pl.Int64, "Scan": pl.Int64, "MSLevel": pl.Int64, "RT": pl.Float64,
        "PrecursorMass": pl.Float64, "#Peaks": pl.Int64,
    }), "peaks": peaks, "ms1": ms1, "chromatogram": chromatogram}


def spectra(peaks: pl.DataFrame) -> pl.DataFrame:
    """One row per spectrum in the columns of the library's raw spectrum
    (``ms_data.spectrum``): its peaks as ``MonoMass_Anno``/``SumIntensity_Anno``
    and, since OpenDDA does not deconvolve, empty mass lists."""
    return (peaks.group_by("index", maintain_order=True)
            .agg(pl.col("mz").alias("MonoMass_Anno"), pl.col("intensity").alias("SumIntensity_Anno"))
            .with_columns(pl.lit([], dtype=pl.List(pl.Float64)).alias("MonoMass"),
                          pl.lit([], dtype=pl.List(pl.Float64)).alias("SumIntensity")))


def state(ctx):
    """The Insight state the components of one Single run share (a picked
    spectrum stays picked across its views)."""
    from src.layout_manager import experiment_state

    return experiment_state("viewer-file", ctx.index, file_id(ctx.run_id, ctx.file))


class Shown:
    """The run (mzML file) a Single run component draws, with where its
    components keep their caches."""

    def __init__(self, ctx, path: Path, whole: bool, results_dir: Path):
        self.ctx, self.path, self.whole = ctx, path, whole
        self.stem = ctx.file
        self.cache_dir = Path(results_dir, "insight_cache", f"file_{self.stem}")
        self.version = int(Path(path).stat().st_mtime)

    def cache_id(self, name: str) -> str:
        return f"{name}_{self.stem}_{self.version}"

    def tables(self) -> dict[str, pl.DataFrame]:
        return _tables(str(self.path), self.version)

    def has_ms1(self) -> bool:
        """Whether the file holds MS1 spectra (only the identified spectra,
        all MS2, are there when the input file is not in the workspace)."""
        return not self.tables()["chromatogram"].is_empty()

    def no_ms1(self, title: str, help_text=None) -> None:
        """Say that a component of the MS1 spectra has nothing to show."""
        import streamlit as st

        st.markdown(f"**{title}**", help=help_text)
        st.info("This run's file holds no MS1 spectra: the input file is not in this workspace, "
                "only the spectra the run identified.")

    def note(self) -> None:
        """Say when only the identified spectra of the file are there."""
        import streamlit as st

        if not self.whole:
            st.caption("The input file is not in this workspace: only the spectra the run identified are shown.")


def shown(ctx) -> Shown | None:
    """What a Single run component of ``ctx`` draws; None (with a note) when
    the file's spectra are nowhere to be found."""
    import streamlit as st

    from src.common.results_helpers import get_results_dir, get_run_params

    if not getattr(ctx, "file", None):
        st.info("Pick a run under Single run.")
        return None
    results_dir = get_results_dir(ctx.workspace, ctx.run_id)
    path, whole = mzml_path(get_run_params(ctx.workspace, ctx.run_id), results_dir, ctx.file)
    if path is None:
        st.info("The spectra of this run are not in this workspace.")
        return None
    return Shown(ctx, path, whole, results_dir)


def _tables(path: str, version: int) -> dict[str, pl.DataFrame]:
    import streamlit as st

    @st.cache_data(show_spinner="Reading the spectra of this run...", max_entries=4)
    def read(path: str, version: int):
        return read_tables(Path(path))

    return read(path, version)
