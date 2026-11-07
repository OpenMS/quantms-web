"""Tests for the run overview (key figures and per-file table)."""
import pandas as pd
import pyopenms as poms
import pytest

# src.common.overview is imported inside the tests, not at collection: importing
# it pulls in src.workflow.ParameterManager with the real streamlit, and
# test_parameter_defaults (collected later) needs to import it with a mock.


def _mzml(path, n_ms1=2, n_ms2=5, n_profile=0):
    exp = poms.MSExperiment()
    for level, n in ((1, n_ms1), (2, n_ms2)):
        for i in range(n):
            spec = poms.MSSpectrum()
            spec.setMSLevel(level)
            if level == 2:
                spec.setType(2 if i < n_profile else 1)
            spec.setRT(float(i))
            spec.set_peaks(([100.0, 200.0], [1.0, 2.0]))
            exp.addSpectrum(spec)
    poms.MzMLFile().store(str(path), exp)


def _idxml(path, hits):
    """hits: (sequence, charge, observed m/z, "target" or "decoy"[, protein accession])."""
    prot = poms.ProteinIdentification()
    prot.setIdentifier("run")
    proteins = []
    for _, _, _, _, *accession in hits:
        if accession:
            protein = poms.ProteinHit()
            protein.setAccession(accession[0])
            proteins.append(protein)
    prot.setHits(proteins)
    peptides = []
    for seq, charge, mz, td, *accession in hits:
        pep = poms.PeptideIdentification()
        pep.setIdentifier("run")
        pep.setMZ(mz)
        hit = poms.PeptideHit()
        hit.setSequence(poms.AASequence.fromString(seq))
        hit.setCharge(charge)
        hit.setMetaValue("target_decoy", td)
        if accession:
            evidence = poms.PeptideEvidence()
            evidence.setProteinAccession(accession[0])
            hit.setPeptideEvidences([evidence])
        pep.setHits([hit])
        peptides.append(pep)
    if hasattr(poms, "PeptideIdentificationList"):
        container = poms.PeptideIdentificationList()
        for pep in peptides:
            container.push_back(pep)
        peptides = container
    path.parent.mkdir(parents=True, exist_ok=True)
    poms.IdXMLFile().store(str(path), [prot], peptides)


def _mz(seq, charge, ppm=0.0, isotope=0):
    theo = poms.AASequence.fromString(seq).getMZ(charge)
    return theo * (1 + ppm / 1e6) + isotope * 1.0033548378 / charge


def test_count_ms2_spectra(tmp_path):
    from src.common.overview import count_ms2_spectra

    _mzml(tmp_path / "a.mzML", n_ms1=3, n_ms2=7)
    assert count_ms2_spectra(tmp_path / "a.mzML") == 7
    assert count_ms2_spectra(tmp_path / "missing.mzML") is None


def test_precursor_error_corrects_isotope_picks():
    from src.common.overview import precursor_error_ppm

    seq = poms.AASequence.fromString("PEPTIDEK")
    assert precursor_error_ppm(seq, 2, _mz("PEPTIDEK", 2, ppm=3.0)) == pytest.approx(3.0, abs=0.01)
    assert precursor_error_ppm(seq, 2, _mz("PEPTIDEK", 2, ppm=-2.0, isotope=1)) == pytest.approx(-2.0, abs=0.05)
    assert precursor_error_ppm(seq, 0, 500.0) is None


def test_search_times_take_the_last_run(tmp_path):
    from src.common.overview import search_times

    log = tmp_path / "commands-and-run-times.log"
    log.write_text(
        "Process finished:\nCometAdapter -in /data/a.mzML -out x.idXML\nTotal time to run command: 10.00 seconds\n\n"
        "Process finished:\nPercolatorAdapter -in x.idXML\nTotal time to run command: 3.00 seconds\n\n"
        "Process finished:\nProSE -in /data/a.mzML -out_idxml x.idXML\nTotal time to run command: 12.50 seconds\n\n"
    )
    assert search_times(log, ["/data/a.mzML", "/data/b.mzML"]) == {"a.mzML": 12.5}
    assert search_times(tmp_path / "none.log", ["/data/a.mzML"]) == {}


def test_run_overview(tmp_path):
    from src.common.overview import run_overview

    wf = tmp_path / "topp-workflow"
    results = wf / "results"
    a, b = tmp_path / "a.mzML", tmp_path / "b.mzML"
    _mzml(a, n_ms2=4)
    _mzml(b, n_ms2=10)
    _idxml(results / "psm_filter" / "a_filter.idXML", [
        ("PEPTIDEK", 2, _mz("PEPTIDEK", 2, ppm=2.0), "target"),
        ("PEPTIDEK", 3, _mz("PEPTIDEK", 3, ppm=4.0), "target"),
        ("ELVISK", 2, _mz("ELVISK", 2, ppm=-1.0), "decoy"),
    ])
    _idxml(results / "percolator_results" / "a_per.idXML", [("PEPTIDEK", 2, 500.0, "target")])
    (wf / "logs").mkdir(parents=True)
    (wf / "logs" / "commands-and-run-times.log").write_text(
        f"Process finished:\nProSE -in {a} -out_idxml x\nTotal time to run command: 5.00 seconds\n\n"
    )

    per_file, figures = run_overview(wf, [str(a), str(b)])

    row_a, row_b = per_file.to_dict("records")
    assert row_a["MS2 spectra"] == 4
    assert row_a["Target PSMs"] == 2 and row_a["Decoy PSMs"] == 1
    assert row_a["Peptides"] == 1
    assert row_a["ID rate (%)"] == 50.0
    assert row_a["Mass error (ppm)"] == pytest.approx(3.0, abs=0.05)
    assert row_a["Percolator"] and row_a["Search (s)"] == 5
    # A file the run did not reach shows empty cells, not zeros
    assert row_b["MS2 spectra"] == 10 and pd.isna(row_b["Target PSMs"])
    assert not row_b["Percolator"]

    assert figures["Target PSMs"] == 2
    assert figures["Peptides"] == 1
    assert figures["PSMs per MS2 (%)"] == round(100 * 2 / 14, 1)
    assert figures["Files"] == 2


def test_run_overview_file_without_psms(tmp_path):
    """A file whose PSMs were all filtered out counts as zero, not as a crash."""
    from src.common.overview import run_overview

    wf = tmp_path / "topp-workflow"
    a, b = tmp_path / "a.mzML", tmp_path / "b.mzML"
    _mzml(a, n_ms2=4)
    _mzml(b, n_ms2=4)
    _idxml(wf / "results" / "psm_filter" / "a_filter.idXML", [("PEPTIDEK", 2, _mz("PEPTIDEK", 2), "target")])
    _idxml(wf / "results" / "psm_filter" / "b_filter.idXML", [])

    per_file, figures = run_overview(wf, [str(a), str(b)])

    row_a, row_b = per_file.to_dict("records")
    assert row_a["Target PSMs"] == 1
    assert row_b["Target PSMs"] == 0 and row_b["Peptides"] == 0
    assert figures["Target PSMs"] == 1


def test_run_overview_diagnosis_columns(tmp_path):
    """Profile spectra, contaminant PSMs (trypsin left out) and PSMs at the window edge."""
    from src.common.overview import run_overview

    wf = tmp_path / "topp-workflow"
    a = tmp_path / "a.mzML"
    _mzml(a, n_ms2=4, n_profile=1)
    _idxml(wf / "results" / "psm_filter" / "a_filter.idXML", [
        ("PEPTIDEK", 2, _mz("PEPTIDEK", 2, ppm=19.0), "target", "CON__sp|P04264|K2C1_HUMAN"),
        ("PEPTIDEK", 3, _mz("PEPTIDEK", 3, ppm=1.0), "target", "CON__sp|P00761|TRYP_PIG"),
        ("ELVISK", 2, _mz("ELVISK", 2, ppm=-2.0), "target", "sp|P12345|ABC_HUMAN"),
        ("ELVISR", 2, _mz("ELVISR", 2, ppm=-1.0), "target"),
    ])

    row = run_overview(wf, [str(a)], precursor_ppm=20.0)[0].iloc[0]
    assert row["Profile MS2 (%)"] == 25.0
    assert row["Contaminant PSMs (%)"] == 25.0
    assert row["At window edge (%)"] == 25.0
    # Without a ppm tolerance the edge stays empty
    assert pd.isna(run_overview(wf, [str(a)])[0].iloc[0]["At window edge (%)"])
