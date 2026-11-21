"""Splitting one Percolator result over all files back into one idXML per file."""

from pathlib import Path

import pytest

oms = pytest.importorskip("pyopenms")

from src.common.pooled_fdr import split_pooled, _load  # noqa: E402


def _peptide_list(peps):
    container = oms.PeptideIdentificationList()
    for pep in peps:
        container.push_back(pep)
    return container


def _psm(identifier, reference, sequence, score):
    pep = oms.PeptideIdentification()
    pep.setIdentifier(identifier)
    pep.setSpectrumReference(reference)
    pep.setMZ(500.0)
    pep.setRT(100.0)
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString(sequence))
    hit.setScore(score)
    hit.setCharge(2)
    pep.setHits([hit])
    return pep


def _search_result(path, identifier, ms_run, psms):
    run = oms.ProteinIdentification()
    run.setIdentifier(identifier)
    run.setSearchEngine("Sage")
    run.setPrimaryMSRunPath([ms_run])
    hit = oms.ProteinHit()
    hit.setAccession("P1")
    run.setHits([hit])
    oms.IdXMLFile().store(str(path), [run], _peptide_list(psms))


def _merged(path, files):
    """What PercolatorAdapter writes for several inputs."""
    run = oms.ProteinIdentification()
    run.setIdentifier("TopPerc_2026")
    run.setSearchEngine("Sage")
    run.setMetaValue("percolator_version", "3.07")
    peps = []
    for index, psms in enumerate(files):
        for reference, sequence, pep_value in psms:
            peps.append(_psm("TopPerc_2026", f"file={index},{reference}", sequence, pep_value))
    oms.IdXMLFile().store(str(path), [run], _peptide_list(peps))


def test_split_restores_files(tmp_path: Path):
    a, b = tmp_path / "a.idXML", tmp_path / "b.idXML"
    _search_result(a, "run_a", "a.mzML", [_psm("run_a", "scan=1", "PEPTIDE", 10.0)])
    _search_result(b, "run_b", "b.mzML", [
        _psm("run_b", "scan=7", "PEPTIDEK", 9.0),
        _psm("run_b", "scan=8", "ELVISK", 8.0),
    ])
    merged = tmp_path / "pooled.idXML"
    _merged(merged, [
        [("scan=1", "PEPTIDE", 0.01)],
        [("scan=7", "PEPTIDEK", 0.02), ("scan=8", "ELVISK", 0.5)],
    ])
    out_a, out_b = tmp_path / "out" / "a_per.idXML", tmp_path / "out" / "b_per.idXML"

    assert split_pooled(str(merged), [str(a), str(b)], [str(out_a), str(out_b)]) == [1, 2]

    proteins, peptides = _load(str(out_b))
    # idXML writes identifiers of its own; the matches must point at the run
    assert proteins[0].getMetaValue("percolator_version") == "3.07"
    paths = proteins[0].getPrimaryMSRunPath()
    assert [p.decode() if isinstance(p, bytes) else p for p in paths] == ["b.mzML"]
    assert [p.getSpectrumReference() for p in peptides] == ["scan=7", "scan=8"]
    assert {p.getIdentifier() for p in peptides} == {proteins[0].getIdentifier()}
    assert [p.getHits()[0].getScore() for p in peptides] == pytest.approx([0.02, 0.5])

    proteins, peptides = _load(str(out_a))
    paths = proteins[0].getPrimaryMSRunPath()
    assert [p.decode() if isinstance(p, bytes) else p for p in paths] == ["a.mzML"]
    assert [p.getSpectrumReference() for p in peptides] == ["scan=1"]


def test_index_fallback_becomes_empty(tmp_path: Path):
    a, b = tmp_path / "a.idXML", tmp_path / "b.idXML"
    _search_result(a, "run_a", "a.mzML", [])
    _search_result(b, "run_b", "b.mzML", [])
    merged = tmp_path / "pooled.idXML"
    _merged(merged, [[("index=1", "PEPTIDE", 0.01)], []])
    out = [str(tmp_path / "a_out.idXML"), str(tmp_path / "b_out.idXML")]

    assert split_pooled(str(merged), [str(a), str(b)], out) == [1, 0]
    _, peptides = _load(out[0])
    assert peptides[0].getSpectrumReference() == ""
    assert Path(out[1]).exists()


def test_match_without_file_position_is_refused(tmp_path: Path):
    a = tmp_path / "a.idXML"
    _search_result(a, "run_a", "a.mzML", [])
    merged = tmp_path / "pooled.idXML"
    run = oms.ProteinIdentification()
    run.setIdentifier("x")
    oms.IdXMLFile().store(str(merged), [run], _peptide_list([_psm("x", "scan=1", "PEPTIDE", 0.1)]))

    with pytest.raises(ValueError, match="file position"):
        split_pooled(str(merged), [str(a)], [str(tmp_path / "out.idXML")])
