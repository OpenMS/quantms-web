"""Explaining a run that ended without PSMs or proteins."""
import pyopenms as poms

from src.common.empty_result import diagnose_empty_run, fdr_summary, run_fdr_levels


def _idxml(path, hits):
    """hits: ("target" or "decoy", q-value) per spectrum."""
    prot = poms.ProteinIdentification()
    prot.setIdentifier("run")
    peptides = []
    for td, q in hits:
        pep = poms.PeptideIdentification()
        pep.setIdentifier("run")
        hit = poms.PeptideHit()
        hit.setSequence(poms.AASequence.fromString("PEPTIDEK"))
        hit.setCharge(2)
        hit.setMetaValue("target_decoy", td)
        if q is not None:
            hit.setMetaValue("MS:1001491", q)
        pep.setHits([hit])
        peptides.append(pep)
    if hasattr(poms, "PeptideIdentificationList"):
        container = poms.PeptideIdentificationList()
        for pep in peptides:
            container.push_back(pep)
        peptides = container
    path.parent.mkdir(parents=True, exist_ok=True)
    poms.IdXMLFile().store(str(path), [prot], peptides)


def _run(tmp_path, comet, percolator=None, filtered=None, log=""):
    run = tmp_path / "run"
    _idxml(run / "results" / "comet_results" / "a_comet.idXML", comet)
    if percolator is not None:
        _idxml(run / "results" / "percolator_results" / "a_per.idXML", percolator)
    if filtered is not None:
        _idxml(run / "results" / "psm_filter" / "a_filter.idXML", filtered)
    (run / "logs").mkdir(parents=True)
    (run / "logs" / "minimal.log").write_text(log)
    return run


LFQ = {"psm-fdr-percent": 1.0, "analysis-mode": "LFQ"}


def test_fdr_levels_and_not_applied():
    assert run_fdr_levels(LFQ) == {"psm": 1.0, "protein": 5.0}  # ProteomicsLFQ's own default
    params = LFQ | {"_defaults": {"ProteomicsLFQ": {"proteinFDR": 0.01}}}
    assert run_fdr_levels(params)["protein"] == 1.0
    params["ProteomicsLFQ"] = {"proteinFDR": 1.0}  # the user switched it off
    assert run_fdr_levels(params)["protein"] is None
    assert fdr_summary(params) == "PSM FDR 1 % (q ≤ 0.01) · protein FDR not applied"
    # Isobaric runs have no protein FDR filter; 100 % switches the PSM filter off
    tmt = {"psm-fdr-percent": 100.0, "analysis-mode": "TMT"}
    assert fdr_summary(tmt) == "PSM FDR not applied · protein FDR not applied"


def test_run_with_psms_is_not_empty(tmp_path):
    run = _run(tmp_path, [("target", None)], [("target", 0.001)], [("target", 0.001)])
    assert diagnose_empty_run(run, LFQ) is None


def test_fdr_filter_left_nothing(tmp_path):
    run = _run(
        tmp_path,
        comet=[("target", None)] * 3 + [("decoy", None)],
        percolator=[("target", 0.04), ("target", 0.08), ("target", 0.2), ("decoy", 0.5)],
        filtered=[("decoy", 0.5)],  # a decoy left over is still no result
    )
    calls = []
    d = diagnose_empty_run(run, LFQ, lambda: calls.append(1) or 400)
    assert d["what"] == "PSM" and d["step"] == "PSM FDR filter"
    assert d["psm_fdr"] == "1 % (q ≤ 0.01)"
    assert [c["targets"] for c in d["counts"]] == [3, 3, 0]
    assert "3 target PSMs reached the filter; the best had q = 0.04" in d["causes"][0]
    assert "Only 400 MS2 spectra" in d["causes"][1]
    assert calls == [1]


def test_search_found_nothing(tmp_path):
    run = _run(tmp_path, comet=[("decoy", None)])
    d = diagnose_empty_run(run, LFQ, 50_000)
    assert d["step"] == "database search"
    assert [c["targets"] for c in d["counts"]] == [0]
    assert "species" in d["causes"][0]


def test_protein_fdr_left_nothing(tmp_path):
    run = _run(
        tmp_path, [("target", None)], [("target", 0.001)], [("target", 0.001)],
        log="Error: Unexpected internal error (No proteins left after FDR filtering. ...)",
    )
    d = diagnose_empty_run(run, LFQ | {"ProteomicsLFQ": {"proteinFDR": 0.01}})
    assert d["what"] == "protein"
    assert d["protein_fdr"] == "1 % (q ≤ 0.01)"
    assert d["counts"][-1] == {"step": "protein FDR filter", "targets": 0, "best_q": None, "unit": "proteins"}


def test_run_that_crashed_mid_way_is_not_explained_as_empty(tmp_path):
    # PSMs after rescoring, then the run stopped before the filter wrote anything
    run = _run(tmp_path, [("target", None)], [("target", 0.001)])
    assert diagnose_empty_run(run, LFQ) is None
