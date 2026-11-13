import pandas as pd
import pyopenms as poms
import pytest

from src.common.charge_cleavage import ALL_FILES, read_psms, searched_enzyme, summary

try:  # pyOpenMS >= 3.5.0
    from pyopenms import PeptideIdentificationList
except ImportError:  # pyOpenMS <= 3.4.x
    PeptideIdentificationList = None


def _psm(sequence, charge, decoy=False):
    hit = poms.PeptideHit()
    hit.setSequence(poms.AASequence.fromString(sequence))
    hit.setCharge(charge)
    hit.setScore(0.01)
    hit.setMetaValue("target_decoy", "decoy" if decoy else "target")
    pep = poms.PeptideIdentification()
    pep.setScoreType("q-value")
    pep.setMZ(500.0)
    pep.setRT(100.0)
    pep.setHits([hit])
    return pep


def _write_idxml(path, peptides, enzyme="Trypsin"):
    prot = poms.ProteinIdentification()
    params = prot.getSearchParameters()
    params.digestion_enzyme = poms.ProteaseDB().getEnzyme(enzyme)
    prot.setSearchParameters(params)
    if PeptideIdentificationList:
        pep_list = PeptideIdentificationList()
        for p in peptides:
            pep_list.push_back(p)
    else:
        pep_list = peptides
    poms.IdXMLFile().store(str(path), [prot], pep_list)


def test_read_psms_counts_missed_cleavages_and_skips_decoys(tmp_path):
    path = tmp_path / "run1_prose_perc_filter.idXML"
    _write_idxml(path, [
        _psm("PEPTIDEK", 2),
        _psm("PEPKTIDEK", 3),  # one missed cleavage
        _psm("PEPKPTIDEK", 2),  # KP is no site for Trypsin
        _psm("PEPRTIDEKAK", 4),  # two
        _psm("KEDITPEP", 2, decoy=True),
    ])
    psms, enzyme = read_psms(path, ["run1.mzML"])
    assert enzyme == "Trypsin"
    assert psms["file"].unique().tolist() == ["run1.mzML"]
    assert psms["missed_cleavages"].tolist() == [0, 1, 0, 2]
    assert psms["charge"].tolist() == [2, 3, 2, 4]


def test_unspecific_search_has_no_missed_cleavages(tmp_path):
    path = tmp_path / "run1_prose_perc_filter.idXML"
    _write_idxml(path, [_psm("PEPKTIDEK", 2)], enzyme="unspecific cleavage")
    psms, enzyme = read_psms(path, ["run1.mzML"])
    assert enzyme == "unspecific cleavage"
    assert psms["missed_cleavages"].isna().all()
    per_file, charges, missed = summary(psms)
    assert missed.empty
    assert "Mean missed" not in per_file.columns
    assert per_file["Charge 2 (%)"].tolist() == [100.0, 100.0]


def test_searched_enzyme_defaults_to_trypsin():
    assert searched_enzyme([]) == "Trypsin/P"
    assert searched_enzyme([poms.ProteinIdentification()]) == "Trypsin/P"


def test_summary_per_file_and_all_files(tmp_path):
    a, b = tmp_path / "a_filter.idXML", tmp_path / "b_filter.idXML"
    _write_idxml(a, [_psm("PEPTIDEK", 2), _psm("PEPTIDEK", 2), _psm("PEPKTIDEK", 3), _psm("PEPTIDER", 1)])
    _write_idxml(b, [_psm("PEPTIDEK", 5), _psm("PEPKTIDEKAR", 2)])
    psms = pd.concat([read_psms(p, ["a.mzML", "b.mzML"])[0] for p in (a, b)], ignore_index=True)
    per_file, charges, missed = summary(psms)
    assert per_file["File"].tolist() == [ALL_FILES, "a.mzML", "b.mzML"]
    assert per_file["file_id"].tolist() == [0, 1, 2]
    assert per_file["Target PSMs"].tolist() == [6, 4, 2]
    row_a = per_file[per_file["File"] == "a.mzML"].iloc[0]
    assert row_a["Charge 1 (%)"] == 25.0 and row_a["Charge 2 (%)"] == 50.0 and row_a["Charge 3 (%)"] == 25.0
    # Charge 5 is counted with 4+
    assert per_file[per_file["File"] == "b.mzML"].iloc[0]["Charge 4+ (%)"] == 50.0
    # Peptides, not PSMs: PEPTIDEK twice in a counts once (3 peptides: 0, 1, 0 missed)
    assert row_a["0 missed (%)"] == pytest.approx(66.7)
    assert row_a["Mean missed"] == pytest.approx(0.33)
    assert per_file[per_file["File"] == "b.mzML"].iloc[0]["2+ missed (%)"] == 50.0
    shares = charges[charges["file_id"] == 1].set_index("bin")["share"]
    assert shares.sum() == pytest.approx(100.0)
    assert sorted(charges["bin"].unique()) == [1, 2, 3, 4]
    assert sorted(missed["bin"].unique()) == [0, 1, 2]
