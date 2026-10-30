import math

import pandas as pd
import pyopenms as poms
import pytest

from src.common.qc_plots import (
    C13_C12_DIFF,
    precursor_error_histogram,
    precursor_error_ppm,
    psm_qc_table,
    psms_per_file,
    score_distribution,
)

try:  # pyOpenMS >= 3.5.0
    from pyopenms import PeptideIdentificationList
except ImportError:  # pyOpenMS <= 3.4.x
    PeptideIdentificationList = None


def _psm(sequence, charge, mz_offset_ppm, decoy, svm, q, isotope=0):
    seq = poms.AASequence.fromString(sequence)
    theo = seq.getMZ(charge)
    hit = poms.PeptideHit()
    hit.setSequence(seq)
    hit.setCharge(charge)
    hit.setScore(0.01)
    hit.setMetaValue("target_decoy", "decoy" if decoy else "target")
    hit.setMetaValue("MS:1001492", svm)
    hit.setMetaValue("MS:1001491", q)
    pep = poms.PeptideIdentification()
    pep.setScoreType("pep")
    pep.setMZ(theo * (1 + mz_offset_ppm * 1e-6) + isotope * C13_C12_DIFF / charge)
    pep.setRT(100.0)
    pep.setHits([hit])
    return pep


def _write_idxml(path, peptides):
    if PeptideIdentificationList:
        pep_list = PeptideIdentificationList()
        for p in peptides:
            pep_list.push_back(p)
    else:
        pep_list = peptides
    poms.IdXMLFile().store(str(path), [poms.ProteinIdentification()], pep_list)


def test_precursor_error_removes_isotope_error():
    theo = 500.0
    assert precursor_error_ppm(theo * (1 + 5e-6), theo, 2) == pytest.approx(5.0, abs=1e-6)
    shifted = theo * (1 + 5e-6) + C13_C12_DIFF / 2
    assert precursor_error_ppm(shifted, theo, 2) == pytest.approx(5.0, abs=1e-6)
    assert math.isnan(precursor_error_ppm(theo, theo, 0))


def test_psm_qc_table_reads_hits(tmp_path):
    path = tmp_path / "run1_comet_perc.idXML"
    _write_idxml(path, [
        _psm("PEPTIDEK", 2, 3.0, False, 2.5, 0.001),
        _psm("PEPTIDER", 2, -4.0, False, 1.0, 0.02, isotope=1),
        _psm("KEDITPEP", 2, 0.0, True, -1.0, 0.5),
    ])
    df = psm_qc_table(path, ["run1.mzML"])
    assert list(df["file"].unique()) == ["run1.mzML"]
    assert df["decoy"].tolist() == [False, False, True]
    assert df["percolator_score"].tolist() == [2.5, 1.0, -1.0]
    assert df["q_value"].tolist() == pytest.approx([0.001, 0.02, 0.5])
    assert df["precursor_error_ppm"].tolist()[:2] == pytest.approx([3.0, -4.0], abs=1e-3)


def test_figures_build(tmp_path):
    df = pd.DataFrame({
        "file": ["a.mzML", "a.mzML", "b.mzML", "b.mzML"],
        "decoy": [False, True, False, False],
        "score": [0.01, 0.9, 0.02, 0.03],
        "score_type": ["pep"] * 4,
        "percolator_score": [2.0, -1.0, 1.5, 0.5],
        "q_value": [0.001, 0.5, 0.005, 0.02],
        "charge": [2, 2, 3, 2],
        "precursor_error_ppm": [1.0, 5.0, -2.0, 0.5],
    })
    fig = score_distribution(df, fdr=0.01)
    assert fig.layout.xaxis.title.text == "Percolator score"
    assert [t.name for t in fig.data] == ["Targets (3)", "Decoys (1)"]
    # the cutoff line sits at the lowest passing target score
    assert fig.layout.shapes[0].x0 == 1.5

    err = precursor_error_histogram(df)
    assert len(err.data[0].x) == 3  # decoys left out

    bars = psms_per_file(df)
    assert list(bars.data[0].y) == ["a.mzML", "b.mzML"]
    assert list(bars.data[0].x) == [1, 2]


def test_score_distribution_falls_back_to_evalue():
    df = pd.DataFrame({
        "file": ["a.mzML"] * 2, "decoy": [False, True], "score": [1e-5, 1.0],
        "score_type": ["expect"] * 2, "percolator_score": [float("nan")] * 2,
        "q_value": [float("nan")] * 2, "charge": [2, 2], "precursor_error_ppm": [0.0, 0.0],
    })
    fig = score_distribution(df)
    assert fig.layout.xaxis.title.text == "-log10(E-value)"
    assert list(fig.data[0].x) == pytest.approx([5.0])
