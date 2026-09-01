from pathlib import Path

import pyopenms as poms

from src.common.results_helpers import (
    extract_filename_from_idxml,
    extract_scan_from_ref,
    extract_scan_number,
    get_workflow_dir,
    load_idxml,
)


def test_get_workflow_dir():
    assert get_workflow_dir("/tmp/ws") == Path("/tmp/ws", "topp-workflow")


def test_extract_scan_from_ref():
    assert extract_scan_from_ref("controllerType=0 controllerNumber=1 scan=1234") == 1234
    assert extract_scan_from_ref("no scan token here") == 0


def test_extract_scan_number():
    assert extract_scan_number("file=abc scan=42") == 42


def test_extract_filename_from_idxml_strips_suffixes():
    assert extract_filename_from_idxml(Path("02COVID_filter.idXML")) == "02COVID.mzML"
    assert extract_filename_from_idxml(Path("sample_comet.idXML")) == "sample.mzML"
    assert extract_filename_from_idxml(Path("run_per.idXML")) == "run.mzML"


def _write_idxml(path):
    """Write a small idXML with one protein and two peptide hits."""
    prot = poms.ProteinIdentification()
    prot.setIdentifier("SEARCH_1")
    hit = poms.ProteinHit()
    hit.setAccession("sp|P12345|TEST")
    prot.setHits([hit])

    peptides = []
    for i, seq in enumerate(("PEPTIDEK", "ELVISLIVESR")):
        pep = poms.PeptideIdentification()
        pep.setIdentifier("SEARCH_1")
        pep.setRT(100.0 + i)
        pep.setMZ(500.0 + i)
        pep.setMetaValue("spectrum_reference", f"scan={i + 1}")
        pep_hit = poms.PeptideHit()
        pep_hit.setSequence(poms.AASequence.fromString(seq))
        pep_hit.setCharge(2)
        evidence = poms.PeptideEvidence()
        evidence.setProteinAccession("sp|P12345|TEST")
        pep_hit.setPeptideEvidences([evidence])
        pep.setHits([pep_hit])
        peptides.append(pep)

    # pyOpenMS >= 3.5.0 wants the dedicated container here too.
    if hasattr(poms, "PeptideIdentificationList"):
        container = poms.PeptideIdentificationList()
        for pep in peptides:
            container.push_back(pep)
        peptides = container

    poms.IdXMLFile().store(str(path), [prot], peptides)


def test_load_idxml_reads_identifications(tmp_path):
    """Regression test for pyOpenMS 3.5.0.

    3.5.0 changed the third parameter of ``IdXMLFile.load()`` from a plain
    ``libcpp_vector[PeptideIdentification]`` to a ``PeptideIdentificationList``,
    so passing ``[]`` raised ``can not handle type of (<path>, [], [])``.
    """
    idxml = tmp_path / "sample_comet.idXML"
    _write_idxml(idxml)

    proteins, peptides = load_idxml(idxml)

    # A plain list either way, so callers can index and len() it.
    assert isinstance(peptides, list)
    assert isinstance(proteins, list)
    assert len(proteins) == 1
    assert len(peptides) == 2

    sequences = [
        hit.getSequence().toString()
        for pep in peptides
        for hit in pep.getHits()
    ]
    assert sequences == ["PEPTIDEK", "ELVISLIVESR"]
    assert peptides[0].getRT() == 100.0
    assert peptides[0].getHits()[0].getCharge() == 2


def test_load_idxml_accepts_str_path(tmp_path):
    """Callers pass both Path and str; load_idxml must handle either."""
    idxml = tmp_path / "run_per.idXML"
    _write_idxml(idxml)

    assert len(load_idxml(str(idxml))[1]) == 2
