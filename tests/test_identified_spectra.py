"""The identified spectra a run keeps (results/spectra) and the viewer's use of them."""
import json
import zipfile
from pathlib import Path

import pyopenms as oms
import pytest

from src.common.identified_spectra import (
    SPECTRA_DIR, has_spectra, identified_refs, spectra_sources, write_identified_spectra, write_spectra_subset,
)


def _write_mzml(path: Path, scans=(2, 3, 4, 5)) -> None:
    """One MS1 and an MS2 per scan, each MS2 with a precursor."""
    exp = oms.MSExperiment()
    ms1 = oms.MSSpectrum()
    ms1.setMSLevel(1)
    ms1.setNativeID("controllerType=0 controllerNumber=1 scan=1")
    ms1.set_peaks(([400.0, 500.0], [10.0, 20.0]))
    exp.addSpectrum(ms1)
    for scan in scans:
        ms2 = oms.MSSpectrum()
        ms2.setMSLevel(2)
        ms2.setRT(10.0 * scan)
        ms2.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
        precursor = oms.Precursor()
        precursor.setMZ(500.0 + scan)
        precursor.setCharge(2)
        ms2.setPrecursors([precursor])
        ms2.set_peaks(([100.0 + scan, 200.0 + scan], [1.0, 2.0]))
        exp.addSpectrum(ms2)
    oms.MzMLFile().store(str(path), exp)


def _write_idxml(path: Path, scans, without_hits=()) -> None:
    """PSMs on ``scans``; ``without_hits`` adds identifications that have no hit."""
    prot = oms.ProteinIdentification()
    prot.setIdentifier("SEARCH_1")
    protein_hit = oms.ProteinHit()
    protein_hit.setAccession("sp|P12345|TEST")
    prot.setHits([protein_hit])
    peptides = []
    for scan in [*scans, *without_hits]:
        pep = oms.PeptideIdentification()
        pep.setIdentifier("SEARCH_1")
        pep.setRT(10.0 * scan)
        pep.setMZ(500.0 + scan)
        pep.setMetaValue("spectrum_reference", f"controllerType=0 controllerNumber=1 scan={scan}")
        if scan in scans:
            hit = oms.PeptideHit()
            hit.setSequence(oms.AASequence.fromString("PEPTIDEK"))
            hit.setCharge(2)
            hit.setScore(0.01)
            evidence = oms.PeptideEvidence()
            evidence.setProteinAccession("sp|P12345|TEST")
            hit.setPeptideEvidences([evidence])
            pep.setHits([hit])
        peptides.append(pep)
    if hasattr(oms, "PeptideIdentificationList"):
        container = oms.PeptideIdentificationList()
        for pep in peptides:
            container.push_back(pep)
        peptides = container
    oms.IdXMLFile().store(str(path), [prot], peptides)


def _load(path: Path) -> oms.MSExperiment:
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(path), exp)
    return exp


def test_identified_refs_skip_identifications_without_hits(tmp_path):
    idxml = tmp_path / "a_filter.idXML"
    _write_idxml(idxml, scans=[3], without_hits=[4])
    refs, scans = identified_refs(idxml)
    assert refs == {"controllerType=0 controllerNumber=1 scan=3"}
    assert scans == {3}


def test_subset_keeps_only_identified_ms2_unchanged(tmp_path):
    mzml, idxml = tmp_path / "a.mzML", tmp_path / "a_filter.idXML"
    _write_mzml(mzml)
    _write_idxml(idxml, scans=[3, 5], without_hits=[4])

    out = tmp_path / "results" / SPECTRA_DIR / "a.mzML"
    assert write_spectra_subset(mzml, [idxml], out) == 2

    exp = _load(out)
    assert [s.getNativeID() for s in exp] == [
        "controllerType=0 controllerNumber=1 scan=3", "controllerType=0 controllerNumber=1 scan=5",
    ]
    assert all(s.getMSLevel() == 2 for s in exp)
    first = exp[0]
    assert first.getPrecursors()[0].getMZ() == pytest.approx(503.0)
    assert first.getPrecursors()[0].getCharge() == 2
    assert first.getRT() == pytest.approx(30.0)
    mz, intensity = first.get_peaks()
    assert list(mz) == pytest.approx([103.0, 203.0])
    assert list(intensity) == pytest.approx([1.0, 2.0])
    assert not out.with_name("a.mzML.part").exists()


def test_write_identified_spectra_one_file_per_input(tmp_path):
    inputs = tmp_path / "in"
    inputs.mkdir()
    a, b = inputs / "a.mzML", inputs / "b.mzML"
    _write_mzml(a)
    _write_mzml(b)
    ida, idb = tmp_path / "a_filter.idXML", tmp_path / "b_filter.idXML"
    _write_idxml(ida, scans=[2])
    _write_idxml(idb, scans=[4, 5])

    results = tmp_path / "results"
    assert write_identified_spectra([a, b], [ida, idb], results) == {"a.mzML": 1, "b.mzML": 2}
    assert has_spectra(results)
    assert sorted(p.name for p in (results / SPECTRA_DIR).iterdir()) == ["a.mzML", "b.mzML"]
    assert [s.getNativeID() for s in _load(results / SPECTRA_DIR / "b.mzML")] == [
        "controllerType=0 controllerNumber=1 scan=4", "controllerType=0 controllerNumber=1 scan=5",
    ]


def test_spectra_sources_prefer_the_input_then_the_subset(tmp_path):
    results = tmp_path / "results"
    (results / SPECTRA_DIR).mkdir(parents=True)
    present = tmp_path / "present.mzML"
    present.write_text("x")
    (results / SPECTRA_DIR / "present.mzML").write_text("x")
    (results / SPECTRA_DIR / "moved.mzML").write_text("x")

    sources = spectra_sources([present, tmp_path / "gone" / "moved.mzML", tmp_path / "lost.mzML"], results)
    assert sources == [present, results / SPECTRA_DIR / "moved.mzML", None]


def test_runs_without_spectra_keep_the_input_paths(tmp_path):
    results = tmp_path / "results"
    results.mkdir()
    assert not has_spectra(results)
    assert spectra_sources([tmp_path / "gone.mzML"], results) == [None]


def test_run_zip_carries_the_spectra_but_not_the_caches(tmp_path):
    from src.workflow.RunManager import RunManager

    workflow_dir = tmp_path / "wf"
    workflow_dir.mkdir()
    (workflow_dir / "params.json").write_text(json.dumps({"mzML-files": ["a.mzML"]}))
    runs = RunManager(workflow_dir)
    run_id = runs.new_run("Run 1", workflow_dir / "params.json")
    results = runs.results_dir(run_id)
    (results / SPECTRA_DIR).mkdir(parents=True)
    (results / SPECTRA_DIR / "a.mzML").write_text("<mzML/>")
    (results / "insight_cache").mkdir()
    (results / "insight_cache" / "x.parquet").write_text("x")
    runs.finish(run_id, "finished")

    with zipfile.ZipFile(runs.zip(run_id, exclude=["insight_cache"])) as zf:
        names = zf.namelist()
    assert "results/spectra/a.mzML" in names
    assert not any("insight_cache" in n for n in names)

    other = RunManager(tmp_path / "other")
    imported = other.import_zip(runs.zip(run_id, exclude=["insight_cache"]))
    assert (other.results_dir(imported) / SPECTRA_DIR / "a.mzML").read_text() == "<mzML/>"


def test_rebuild_views_from_the_subset_when_the_input_is_gone(tmp_path):
    pytest.importorskip("openms_insight")
    import polars as pl

    from src.common.psm_views import rebuild_psm_views

    inputs = tmp_path / "in"
    inputs.mkdir()
    mzml = inputs / "a.mzML"
    _write_mzml(mzml)
    results = tmp_path / "results"
    idxml = results / "psm_filter" / "a_filter.idXML"
    idxml.parent.mkdir(parents=True)
    _write_idxml(idxml, scans=[3, 5])
    write_identified_spectra([mzml], [idxml], results)
    mzml.rename(tmp_path / "moved.mzML")

    assert rebuild_psm_views(results, idxml, [mzml], 0.02, False)
    cache = results / "insight_cache"
    assert (cache / "table_a_filter").is_dir()
    peaks = pl.read_parquet(cache / "lineplot_a_filter" / "preprocessed" / "data.parquet")
    assert sorted(set(peaks["scan_id"].to_list())) == [3, 5]
    assert set(peaks["file_index"].to_list()) == {0}
    # Present caches are left alone
    assert rebuild_psm_views(results, idxml, [mzml], 0.02, False)


def test_no_rebuild_for_runs_without_spectra(tmp_path):
    pytest.importorskip("openms_insight")
    from src.common.psm_views import rebuild_psm_views

    mzml = tmp_path / "a.mzML"
    _write_mzml(mzml)
    results = tmp_path / "results"
    idxml = results / "psm_filter" / "a_filter.idXML"
    idxml.parent.mkdir(parents=True)
    _write_idxml(idxml, scans=[3])

    assert not rebuild_psm_views(results, idxml, [mzml], 0.02, False)
    assert not (results / "insight_cache").exists()


def test_search_hits_sort_best_first():
    from src.common.psm_views import score_sort

    # ProSE's ln(hyperscore): higher is better; Percolator, IDFilter and Comet: lower
    assert score_sort(Path("results/search_results/a_prose.idXML")) == "desc"
    assert score_sort(Path("results/percolator_results/a_prose_perc.idXML")) == "asc"
    assert score_sort(Path("results/psm_filter/a_prose_perc_filter.idXML")) == "asc"
    assert score_sort(Path("results/comet_results/a_comet.idXML")) == "asc"
