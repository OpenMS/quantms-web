"""Predicted spectra for the PSM mirror plot."""
import numpy as np
import pandas as pd
import polars as pl
import pytest

from src.common.ms2_prediction import match_peaks, mirror_frame, spectral_angle

pytest.importorskip("alphabase")


def test_to_peptdeep_maps_openms_modifications():
    from src.common.ms2_prediction import to_peptdeep

    assert to_peptdeep("PEPTIDEK") == ("PEPTIDEK", "", "")
    assert to_peptdeep(".(Acetyl)M(Oxidation)PEC(Carbamidomethyl)K") == (
        "MPECK", "Acetyl@Any_N-term;Oxidation@M;Carbamidomethyl@C", "0;1;4",
    )
    assert to_peptdeep(".(TMT6plex)PEK(TMT6plex)") == ("PEK", "TMT6plex@Any_N-term;TMT6plex@K", "0;3")
    # residue-specific N-terminal modification
    assert to_peptdeep("Q(Gln->pyro-Glu)PEK") == ("QPEK", "Gln->pyro-Glu@Q^Any_N-term", "0")


def test_match_peaks_takes_most_intense_peak_in_tolerance():
    measured_mz = np.array([100.0, 100.01, 200.0, 300.5])
    measured_int = np.array([5.0, 9.0, 3.0, 7.0])
    matched = match_peaks(measured_mz, measured_int, np.array([100.005, 200.03, 300.5]), 0.02, False)
    assert matched.tolist() == [9.0, 0.0, 7.0]
    # 20 ppm of 300.5 is 0.006 Da
    assert match_peaks(measured_mz, measured_int, np.array([300.51]), 20, True).tolist() == [0.0]


def test_spectral_angle_bounds():
    assert spectral_angle(np.array([1.0, 2.0]), np.array([2.0, 4.0])) == pytest.approx(1.0)
    assert spectral_angle(np.array([1.0, 0.0]), np.array([0.0, 1.0])) == pytest.approx(0.0)
    assert spectral_angle(np.zeros(2), np.array([0.0, 1.0])) == 0.0


def test_mirror_frame_keeps_measured_peak_ids_and_scales_both_sides():
    measured = pl.DataFrame({"peak_id": [10, 11], "mass": [100.0, 200.0], "intensity": [50.0, 200.0]})
    predicted = pd.DataFrame({"mz": [150.0], "intensity": [0.5], "ion": ["y1¹⁺"]})
    frame = mirror_frame(measured, predicted)
    top, bottom = frame.filter(pl.col("side") == 0), frame.filter(pl.col("side") == 1)
    assert top["peak_id"].to_list() == [10, 11]
    assert top["intensity"].to_list() == [25.0, 100.0]
    assert bottom["peak_id"].to_list() == [-1]
    assert bottom["intensity"].to_list() == [50.0]
    assert bottom["ion"].to_list() == ["y1¹⁺"]


def test_predict_fragments_labels_b_and_y_ions():
    pytest.importorskip("peptdeep")
    from src.common.ms2_prediction import load_model, predict_fragments

    predicted = predict_fragments(load_model(), "PEPTIDEK", 2)
    assert predicted["intensity"].max() == pytest.approx(1.0)
    ions = dict(zip(predicted["ion"], predicted["mz"]))
    # y1 of a C-terminal lysine; doubly charged fragments need a precursor above charge 2
    assert ions["y1¹⁺"] == pytest.approx(147.1128, abs=1e-3)
    assert not any(ion.endswith("²⁺") for ion in ions)
    assert predict_fragments(load_model(), "PEPTIDEXK", 2) is None
