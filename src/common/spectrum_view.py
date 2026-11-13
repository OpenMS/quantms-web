"""Annotated PSM spectrum, mirrored against the AlphaPeptDeep prediction when possible."""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

import pandas as pd
import polars as pl
import pyarrow.parquet as pq
import pyopenms as poms
import streamlit as st
from openms_insight import LinePlot, MirrorPlot

from src.common.ms2_prediction import (
    cosine_similarity, instrument_of, match_peaks, mirror_frame, peptdeep_available,
    spectral_angle,
)
from src.common.identified_spectra import spectra_sources
from src.common.prediction_worker import PredictionWorker
from src.common.psm_features import features_by_id, features_frame
from src.common.results_helpers import parse_idxml, search_tolerances
from src.common.plot_style import spectrum_styling
from src.viewer.scope import zoom


def fragment_tolerance(workflow_params: dict) -> tuple[float, bool]:
    """The search's fragment tolerance: user value, else the workflow's default for it.

    ``workflow_params`` are the parameters of the run shown (``get_run_params``).
    """
    tolerances = search_tolerances(workflow_params)
    return tolerances["fragment"], tolerances["fragment_ppm"]


@st.cache_resource(show_spinner=False)
def _worker() -> PredictionWorker:
    # The model runs in a process of its own: torch must not be loaded in the app server.
    return PredictionWorker()


@st.cache_data(show_spinner=False)
def _instrument(mzml_path: str) -> str:
    return instrument_of(mzml_path)


@st.cache_data(show_spinner=False)
def _measured(parquet: str, mtime: float, file_index: int, scan: int) -> pl.DataFrame:
    # mtime only keys the cache, so a re-run of the workflow invalidates it.
    # Read with pyarrow: on Windows, polars' parquet scan crashed the app server
    # (access violation in a polars worker thread) once the prediction had run.
    table = pq.read_table(
        parquet, columns=["peak_id", "mass", "intensity"],
        filters=[("file_index", "==", file_index), ("scan_id", "==", scan)],
    )
    return pl.from_arrow(table)


@st.cache_data(show_spinner="Reading identifications...")
def psm_lookup(idxml: str, mtime: float) -> pd.DataFrame:
    """Peptide, charge, precursor m/z and retention time of every hit in an idXML, by the id the run's viewers use."""
    id_df, _ = parse_idxml(Path(idxml))
    return id_df.select("id_idx", "sequence", "charge", "mz", "rt").to_pandas().set_index("id_idx")


@st.cache_data(show_spinner="Predicting the spectrum...")
def _predict(peptide: str, charge: int, instrument: str) -> pd.DataFrame | None:
    return _worker().predict(peptide, charge, instrument)


def svg_name(peptide: str, charge: int, scan: int) -> str:
    """File name (without extension) of a PSM's spectrum download, e.g. ``PEPTIDEK_2+_scan1234``.

    Modifications such as ``M(Oxidation)`` become ``M-Oxidation``; anything
    else a file system might reject is dropped.
    """
    name = re.sub(r"[\[\]()]+", "-", str(peptide)).strip("-.")
    name = re.sub(r"[^A-Za-z0-9_+.-]+", "", name).replace("-.", ".")[:80] or "spectrum"
    return f"{name}_{int(charge)}+_scan{int(scan)}"


@st.cache_data(show_spinner=False)
def modified_residues(peptide: str) -> dict[str, float]:
    """Monoisotopic residue masses of the peptide's modified residues, e.g. ``{"M(Oxidation)": 147.0354}``.

    The spectrum's measure tool matches the m/z difference of two clicked peaks
    to residue masses; with these it also recognises the modified residues of
    the PSM shown. Terminal modifications (no residue of their own) and
    sequences pyOpenMS cannot read give nothing.
    """
    try:
        sequence = poms.AASequence.fromString(str(peptide))
    except Exception:  # a sequence pyOpenMS cannot parse: no extra masses
        return {}
    masses: dict[str, float] = {}
    for i in range(sequence.size()):
        residue = sequence.getResidue(i)
        if not residue.isModified():
            continue
        name = f"{residue.getOneLetterCode()}({residue.getModificationName()})"
        masses[name] = round(residue.getMonoWeight(poms.Residue.ResidueType.Internal), 5)
    return masses


def measure_options(peptide: str, tolerance: float, ppm: bool) -> dict:
    """Arguments of the spectrum plots' measure tool for a PSM: its modified residues and,
    for a tolerance in Da, the search's fragment tolerance (a ppm tolerance keeps the plot's default)."""
    options: dict = {"measure_residues": modified_residues(peptide)}
    if not ppm and tolerance > 0:
        options["measure_tolerance"] = float(tolerance)
    return options


def show_spectrum(cache_dir: Path, run: str, state_manager, *, peptide: str, charge: int, file_index: int,
                  scan: int, seqview_key: str, key: str, tolerance: float, ppm: bool,
                  mzml_path: Path | None = None, height: int = 450) -> None:
    """The measured spectrum of a PSM, annotated by the SequenceView ``seqview_key``.

    Mirrored against the spectrum AlphaPeptDeep predicts for the peptide; the
    plain annotated spectrum when the prediction is unavailable.
    """
    parquet = cache_dir / f"lineplot_{run}" / "preprocessed" / "data.parquet"
    predicted, reason = None, None
    if not peptdeep_available():
        reason = "the prediction software (AlphaPeptDeep) is not installed"
    elif not parquet.is_file():
        reason = "this run has no spectrum data; re-run the analysis to create it"
    else:
        instrument = _instrument(str(mzml_path)) if mzml_path and Path(mzml_path).is_file() else "Lumos"
        try:
            predicted = _predict(peptide, int(charge), instrument)
        except Exception as error:  # a model download or prediction failure must not break the page
            reason = f"the prediction failed ({error})"
        if predicted is None and reason is None:
            reason = "the peptide carries a modification or residue the model does not know"

    measure = measure_options(peptide, tolerance, ppm)
    if predicted is None:
        LinePlot(cache_id=f"lineplot_{run}", cache_path=str(cache_dir), **measure, **zoom(LinePlot, key))(
            key=f"{key}_lineplot", state_manager=state_manager, height=height, sequence_view_key=seqview_key,
        )
        return

    measured = _measured(str(parquet), parquet.stat().st_mtime, int(file_index), int(scan))
    if measured.is_empty():
        st.info(
            "The measured spectrum of this PSM was not found. A run opened from a ZIP without its mzML files "
            "only holds the spectra of its final, filtered PSMs; otherwise re-run the analysis to create it."
        )
        return

    observed = match_peaks(
        measured["mass"].to_numpy(), measured["intensity"].to_numpy(), predicted["mz"].to_numpy(), tolerance, ppm,
    )
    similarity = spectral_angle(observed, predicted["intensity"].to_numpy())
    cosine = cosine_similarity(observed, predicted["intensity"].to_numpy())
    found = int((observed > 0).sum())
    c1, c2, c3 = st.columns(3)
    c1.metric(
        "Similarity to prediction", f"{similarity:.2f}",
        help="Normalised spectral angle between the measured and predicted fragment intensities: 1 means identical, "
             "0 nothing in common. Correct matches usually score above 0.7.",
    )
    c2.metric(
        "Cosine similarity", f"{cosine:.2f}",
        help="Cosine of the measured and predicted fragment intensities, the other common measure: 1 means "
             "identical, 0 nothing in common. It runs higher than the spectral angle for the same pair.",
    )
    c3.metric("Predicted ions found", f"{found} of {len(predicted)}",
              help=f"Predicted fragment ions with a measured peak within ±{tolerance:g} {'ppm' if ppm else 'Da'}.")

    data = mirror_frame(measured, predicted)
    digest = hashlib.md5(
        f"{peptide}|{charge}|{file_index}|{scan}|{parquet.stat().st_mtime}|{instrument}".encode()
    ).hexdigest()[:12]
    state_manager.set_selection("mirror_measured", 0)
    state_manager.set_selection("mirror_predicted", 1)
    MirrorPlot(
        cache_id=f"mirror_{digest}",
        data=data.lazy(),
        cache_path=str(cache_dir / "predicted_spectra"),
        filters_top={"mirror_measured": "side"},
        filter_defaults_top={"mirror_measured": 0},
        filters_bottom={"mirror_predicted": "side"},
        filter_defaults_bottom={"mirror_predicted": 1},
        interactivity={"peak": "peak_id"},
        x_column="mass",
        y_column="intensity",
        highlight_column="highlight",
        annotation_column="ion",
        title_top="Measured",
        title_bottom=f"Predicted by AlphaPeptDeep ({instrument})",
        x_label="m/z",
        y_label="Relative intensity (%)",
        styling=spectrum_styling(),
        svg_filename=svg_name(peptide, charge, scan),
        **measure,
        **zoom(MirrorPlot, f"{key}_mirror"),
    )(key=f"{key}_mirror", state_manager=state_manager, height=height, sequence_view_top_key=seqview_key)


def show_selected_psm_spectrum(cache_dir: Path, idxml: Path, state_manager, *, seqview_key: str,
                               workflow_params: dict, mzml_files: list) -> None:
    """``show_spectrum`` for the PSM selected in a run's PSM table (state: identification, file, spectrum)."""
    psm_id = state_manager.get_selection("identification")
    file_index = state_manager.get_selection("file")
    scan = state_manager.get_selection("spectrum")
    hits = psm_lookup(str(idxml), idxml.stat().st_mtime)
    if psm_id is None or file_index is None or scan is None or int(psm_id) not in hits.index:
        st.info("Select a PSM in the table to see its spectrum next to the predicted one.")
        return
    hit = hits.loc[int(psm_id)]
    tolerance, ppm = fragment_tolerance(workflow_params)
    mzml = mzml_files[int(file_index)] if int(file_index) < len(mzml_files) else None
    show_spectrum(
        cache_dir, idxml.stem, state_manager, peptide=hit["sequence"], charge=int(hit["charge"]),
        file_index=int(file_index), scan=int(scan), seqview_key=seqview_key, key=seqview_key,
        tolerance=tolerance, ppm=ppm, mzml_path=mzml,
    )
    show_psm_features(idxml, int(psm_id))


@st.cache_data(show_spinner="Reading scores and features...")
def _features(idxml: str, mtime: float) -> dict:
    # mtime only keys the cache, so a re-run of the workflow invalidates it
    return features_by_id(Path(idxml))


def show_psm_features(idxml: Path, psm_id: int) -> None:
    """Every score and feature the search and rescoring gave one PSM, in an expander."""
    rows = _features(str(idxml), idxml.stat().st_mtime).get(psm_id)
    if not rows:
        return
    with st.expander(f"All scores and features of this PSM ({len(rows)})"):
        st.caption(
            "What ProSE and Percolator recorded for this match. Percolator combines the ProSE values into its "
            "score; the q-value and PEP come from that score."
        )
        st.dataframe(
            features_frame(rows), hide_index=True, use_container_width=True,
            height=min(38 + 35 * len(rows), 600),
            column_config={
                "Source": st.column_config.TextColumn(width="small"),
                "Feature": st.column_config.TextColumn(width="medium"),
                "Meaning": st.column_config.TextColumn(width="large"),
            },
        )


def psm_mzml_files(workflow_params: dict, results_dir: Path | None = None) -> list[Path | None]:
    """A run's mzML inputs, in the order that numbers them in the viewer caches.

    With ``results_dir``, an input that is not there is replaced by the run's
    own subset of its identified spectra (``results/spectra``), or ``None``
    when the run has none.
    """
    from src.WorkflowTest import WorkflowTest

    wf = WorkflowTest()
    names = workflow_params.get("mzML-files") or []
    files = [Path(f) for f in wf.file_manager.get_files(names)] if names else []
    return spectra_sources(files, results_dir) if results_dir is not None else files


def ensure_psm_views(results_dir: Path, idxml: Path, workflow_params: dict) -> bool:
    """Whether the viewer caches of ``idxml`` exist, making them first for a run opened from a ZIP.

    A run's ZIP leaves the caches out but carries the identified spectra; the
    caches are built again from those (or from the input mzML when it is
    there). Runs without the spectra are left as they are.
    """
    from src.common.psm_views import rebuild_psm_views

    if (Path(results_dir, "insight_cache") / f"table_{Path(idxml).stem}").is_dir():
        return True
    tolerance, ppm = fragment_tolerance(workflow_params)
    with st.spinner("Preparing the spectrum views of this run (first time only)..."):
        return rebuild_psm_views(results_dir, idxml, psm_mzml_files(workflow_params), tolerance, ppm)
