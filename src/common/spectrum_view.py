"""Annotated PSM spectrum, mirrored against the AlphaPeptDeep prediction when possible."""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import polars as pl
import streamlit as st
from openms_insight import LinePlot, MirrorPlot

from src.common.ms2_prediction import (
    instrument_of, load_model, match_peaks, mirror_frame, peptdeep_available, predict_fragments, spectral_angle,
)
from src.common.results_helpers import parse_idxml


def fragment_tolerance(workflow_params: dict) -> tuple[float, bool]:
    """The search's fragment tolerance: user value, else the workflow's default for it.

    ``workflow_params`` are the parameters of the run shown (``get_run_params``).
    """
    comet = "CometAdapter-TMT" if workflow_params.get("analysis-mode", "LFQ") == "TMT" else "CometAdapter"
    comet_params = {**workflow_params.get("_defaults", {}).get(comet, {}), **workflow_params.get(comet, {})}
    return float(comet_params.get("fragment_mass_tolerance", 0.02)), comet_params.get("fragment_error_units", "Da") != "Da"


@st.cache_resource(show_spinner="Loading the spectrum prediction model (first time only)...")
def _model():
    return load_model()


@st.cache_data(show_spinner=False)
def _instrument(mzml_path: str) -> str:
    return instrument_of(mzml_path)


@st.cache_data(show_spinner=False)
def _measured(parquet: str, mtime: float, file_index: int, scan: int) -> pl.DataFrame:
    # mtime only keys the cache, so a re-run of the workflow invalidates it
    return (
        pl.scan_parquet(parquet)
        .filter((pl.col("file_index") == file_index) & (pl.col("scan_id") == scan))
        .select("peak_id", "mass", "intensity")
        .collect()
    )


@st.cache_data(show_spinner="Reading identifications...")
def psm_lookup(idxml: str, mtime: float) -> pd.DataFrame:
    """Peptide and charge of every hit in an idXML, by the id the run's viewers use."""
    id_df, _ = parse_idxml(Path(idxml))
    return id_df.select("id_idx", "sequence", "charge").to_pandas().set_index("id_idx")


@st.cache_data(show_spinner="Predicting the spectrum...")
def _predict(peptide: str, charge: int, instrument: str) -> pd.DataFrame | None:
    return predict_fragments(_model(), peptide, charge, instrument)


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

    if predicted is None:
        st.caption(f"No predicted spectrum for comparison: {reason}.")
        LinePlot(cache_id=f"lineplot_{run}", cache_path=str(cache_dir))(
            key=f"{key}_lineplot", state_manager=state_manager, height=height, sequence_view_key=seqview_key,
        )
        return

    measured = _measured(str(parquet), parquet.stat().st_mtime, int(file_index), int(scan))
    if measured.is_empty():
        st.info("The measured spectrum of this PSM was not found. Re-run the analysis to create it.")
        return

    observed = match_peaks(
        measured["mass"].to_numpy(), measured["intensity"].to_numpy(), predicted["mz"].to_numpy(), tolerance, ppm,
    )
    similarity = spectral_angle(observed, predicted["intensity"].to_numpy())
    found = int((observed > 0).sum())
    c1, c2 = st.columns(2)
    c1.metric(
        "Similarity to prediction", f"{similarity:.2f}",
        help="Normalised spectral angle between the measured and predicted fragment intensities: 1 means identical, "
             "0 nothing in common. Correct matches usually score above 0.7.",
    )
    c2.metric("Predicted ions found", f"{found} of {len(predicted)}",
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
        styling={"unhighlightedColor": "#CCCCCC", "highlightColor": "#E74C3C", "selectedColor": "#F3A712"},
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
    mzml = Path(mzml_files[int(file_index)]) if int(file_index) < len(mzml_files) else None
    show_spectrum(
        cache_dir, idxml.stem, state_manager, peptide=hit["sequence"], charge=int(hit["charge"]),
        file_index=int(file_index), scan=int(scan), seqview_key=seqview_key, key=seqview_key,
        tolerance=tolerance, ppm=ppm, mzml_path=mzml,
    )


def psm_mzml_files(workflow_params: dict) -> list[Path]:
    """A run's mzML inputs, in the order that numbers them in the viewer caches."""
    from src.WorkflowTest import WorkflowTest

    wf = WorkflowTest()
    names = workflow_params.get("mzML-files") or []
    return [Path(f) for f in wf.file_manager.get_files(names)] if names else []
