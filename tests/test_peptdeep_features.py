"""Predicted-spectrum features for Percolator (src/common/peptdeep_features.py).

The model is replaced by a stub that predicts fixed b/y ions, so these tests
need pyopenms but not peptdeep or its downloaded models.
"""
import numpy as np
import pandas as pd
import pytest

pyopenms = pytest.importorskip("pyopenms")

from src.common import peptdeep_features as pf  # noqa: E402


def test_cosine_bounds():
    assert pf.cosine(np.array([1.0, 2.0]), np.array([2.0, 4.0])) == pytest.approx(1.0)
    assert pf.cosine(np.array([1.0, 0.0]), np.array([0.0, 1.0])) == 0.0
    assert pf.cosine(np.zeros(3), np.ones(3)) == 0.0


def test_calibration_set_takes_confident_targets_only():
    psms = pd.DataFrame({
        "predictable": [True] * 40 + [False, False],
        "decoy": [False] * 39 + [True, False, False],
        "expect": [1e-4] * 39 + [1e-6, 1e-6, 1e-6],
    })
    chosen = pf.calibration_set(psms)
    assert len(chosen) == 39
    assert not chosen["decoy"].any() and chosen["predictable"].all()


def test_rt_calibration_ignores_wrong_matches():
    rng = np.random.default_rng(1)
    predicted = pd.Series(rng.uniform(0, 1, 200))
    observed = 3600 * predicted + 300 + rng.normal(0, 5, 200)
    observed.iloc[:20] = rng.uniform(300, 3900, 20)  # wrong identifications
    to_seconds = pf.calibrate_rt(predicted, observed, predicted.index)
    assert to_seconds([0.0, 1.0]) == pytest.approx([300, 3900], abs=30)


def test_fragment_vectors_drop_doubly_charged_ions_of_doubly_charged_precursors():
    result = {
        "precursor_df": pd.DataFrame({"frag_start_idx": [0], "frag_stop_idx": [2]}),
        "fragment_mz_df": pd.DataFrame({t: [100.0, 200.0] for t in pf.FRAG_TYPES}),
        "fragment_intensity_df": pd.DataFrame({t: [0.5, 1.0] for t in pf.FRAG_TYPES}),
    }
    assert len(pf.fragment_vectors(result, 0, 2)[0]) == 4  # b and y, 1+ only
    assert len(pf.fragment_vectors(result, 0, 3)[0]) == 8


class StubModel:
    """Predicts y1..y3 of every peptide at unit intensity and RT = position / 10."""

    def predict_all(self, precursors, predict_items, frag_types, multiprocessing):
        n = len(precursors)
        precursors = precursors.copy()
        precursors["frag_start_idx"] = np.arange(n) * 3
        precursors["frag_stop_idx"] = precursors["frag_start_idx"] + 3
        precursors["rt_pred"] = np.arange(n) / 10
        ions = np.tile([300.0, 400.0, 500.0], n)
        return {
            "precursor_df": precursors,
            "fragment_mz_df": pd.DataFrame({t: ions if t == "y_z1" else 0.0 for t in frag_types}),
            "fragment_intensity_df": pd.DataFrame({t: np.ones(3 * n) if t == "y_z1" else 0.0 for t in frag_types}),
        }


def _write_inputs(tmp_path):
    exp = pyopenms.MSExperiment()
    for i, peaks in enumerate([[300.0, 400.0, 500.0], [310.0, 410.0, 510.0]], start=1):
        spectrum = pyopenms.MSSpectrum()
        spectrum.setMSLevel(2)
        spectrum.setNativeID(f"scan={i}")
        spectrum.setRT(60.0 * i)
        spectrum.set_peaks((peaks, [100.0, 100.0, 100.0]))
        exp.addSpectrum(spectrum)
    mzml = tmp_path / "run.mzML"
    pyopenms.MzMLFile().store(str(mzml), exp)

    protein = pyopenms.ProteinIdentification()
    protein.setIdentifier("search")
    protein.setSearchEngine("Comet")
    protein.setHits([pyopenms.ProteinHit()])
    hits = protein.getHits()
    hits[0].setAccession("P1")
    protein.setHits(hits)
    peptides = []
    for i, (seq, td) in enumerate([("PEPTIDEK", "target"), ("KEDITPEP", "decoy")], start=1):
        hit = pyopenms.PeptideHit()
        hit.setSequence(pyopenms.AASequence.fromString(seq))
        hit.setCharge(2)
        hit.setMetaValue("target_decoy", td)
        hit.setMetaValue("expect", 1e-5)
        evidence = pyopenms.PeptideEvidence()
        evidence.setProteinAccession("P1")
        hit.setPeptideEvidences([evidence])
        pep = pyopenms.PeptideIdentification()
        pep.setIdentifier("search")
        pep.setRT(60.0 * i)
        pep.setMetaValue("spectrum_reference", f"scan={i}")
        pep.setHits([hit])
        peptides.append(pep)
    try:
        container = pyopenms.PeptideIdentificationList()
        for pep in peptides:
            container.push_back(pep)
    except AttributeError:
        container = peptides
    idxml = tmp_path / "run_comet.idXML"
    pyopenms.IdXMLFile().store(str(idxml), [protein], container)
    return idxml, mzml


def test_add_features_annotates_hits_and_registers_them_for_percolator(tmp_path):
    pytest.importorskip("alphabase")  # peptide notation conversion; installed with peptdeep
    from src.common.results_helpers import load_idxml

    idxml, mzml = _write_inputs(tmp_path)
    out = tmp_path / "run_pd.idXML"
    summary = pf.add_features(idxml, mzml, out, StubModel(), tolerance=0.02, ppm=False,
                              instrument="QE", log=lambda _: None)

    assert summary["psms"] == summary["predicted"] == 2
    proteins, peptides = load_idxml(out)
    extra = proteins[0].getSearchParameters().getMetaValue("extra_features")
    assert (extra.decode() if isinstance(extra, bytes) else extra).split(",") == pf.FEATURES
    angles = {
        p.getMetaValue("spectrum_reference"): p.getHits()[0].getMetaValue("peptdeep:spectral_angle")
        for p in peptides
    }
    angles = {(k.decode() if isinstance(k, bytes) else k): v for k, v in angles.items()}
    assert angles["scan=1"] == pytest.approx(1.0)  # measured peaks are the predicted ones
    assert angles["scan=2"] == 0.0  # 10 Da off: nothing matches
    for p in peptides:
        assert p.getHits()[0].metaValueExists("peptdeep:delta_rt")
