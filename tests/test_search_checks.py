"""Checks of the search settings against the runs and the database."""

from src.common import instrument
from src.common.search_checks import (
    CID, ETD, HCD, RunDetails, check, read_activation, read_run, run_summary, summarize_database,
)

ORBITRAP_HCD = {"Search:fragment:mass_tolerance": 20.0, "Search:fragment:mass_tolerance_unit": "ppm"}
ION_TRAP_CID = {"Search:fragment:mass_tolerance": 0.5, "Search:fragment:mass_tolerance_unit": "Da"}


def _mzml(activation: str, analyzer: str = "orbitrap") -> str:
    spectrum = (
        '<spectrum><precursorList><precursor><activation>'
        f'<cvParam cvRef="MS" accession="MS:1" name="{activation}" />'
        '<cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="30"/>'
        '</activation></precursor></precursorList></spectrum>'
    )
    return (
        '<mzML><instrumentConfigurationList count="1"><instrumentConfiguration id="IC1">'
        '<cvParam cvRef="MS" accession="MS:1002523" name="Q Exactive HF"/>'
        f'<componentList><analyzer order="2"><cvParam cvRef="MS" accession="MS:2" name="{analyzer}"/>'
        '</analyzer></componentList></instrumentConfiguration></instrumentConfigurationList>'
        + spectrum * 3 + '</mzML>'
    )


def test_activation_is_read_from_the_spectra():
    assert read_activation(_mzml("beam-type collision-induced dissociation")) == HCD
    assert read_activation(_mzml("collision-induced dissociation")) == CID
    assert read_activation(_mzml("electron transfer dissociation")) == ETD
    assert read_activation('<cvParam name="filter string" value="FTMS + p NSI d Full ms2 500.0@etd30.00"/>') == ETD
    assert read_activation("<mzML/>") is None


def test_a_run_file_gives_model_analyzer_and_activation(tmp_path):
    path = tmp_path / "a.mzML"
    path.write_text(_mzml("beam-type collision-induced dissociation"))
    run = read_run(path)
    assert (run.model, run.analyzer, run.activation) == ("Q Exactive HF", instrument.ORBITRAP, HCD)
    assert read_run(tmp_path / "b.raw", tmp_path) is None  # not converted yet
    assert run_summary([run, run]) == "2 × Q Exactive HF (HCD, MS2 in the Orbitrap)"


def test_ion_trap_ms2_with_a_ppm_fragment_tolerance_is_flagged():
    runs = [RunDetails("a", analyzer=instrument.ION_TRAP, activation=CID)]
    notices = check(runs, ORBITRAP_HCD, None)
    assert [n.id for n in notices] == ["fragment-ion-trap"]
    assert notices[0].preset == instrument.ION_TRAP
    assert check(runs, ION_TRAP_CID, None) == []


def test_high_resolution_ms2_with_a_wide_da_tolerance_is_flagged():
    runs = [RunDetails("a", analyzer=instrument.ORBITRAP, activation=HCD)]
    notices = check(runs, ION_TRAP_CID, None)
    assert [(n.id, n.preset) for n in notices] == [("fragment-high-res", instrument.ORBITRAP)]
    assert check(runs, ORBITRAP_HCD, None) == []


def test_etd_and_mixed_fragmentation():
    etd = check([RunDetails("a", activation=ETD), RunDetails("b", activation=HCD)], ORBITRAP_HCD, None)
    assert [n.id for n in etd] == ["fragmentation-etd"]
    mixed = check([RunDetails("a", activation=CID), RunDetails("b", activation=HCD)], ORBITRAP_HCD, None)
    assert [n.id for n in mixed] == ["fragmentation-mixed"]


def _fasta(path, entries):
    path.write_text("".join(f">{h}\nPEPTIDEK\n" for h in entries))
    return path


def test_database_summary_counts_targets_decoys_contaminants_and_organisms(tmp_path):
    path = _fasta(tmp_path / "db.fasta", [
        "sp|P1|A_YEAST Protein A OS=Saccharomyces cerevisiae OX=559292 GN=A",
        "sp|P2|B_YEAST Protein B OS=Saccharomyces cerevisiae OX=559292",
        "sp|P3|C_HUMAN Protein C OS=Homo sapiens OX=9606",
        "CON__P4 keratin",
        "DECOY_sp|P1|A_YEAST",
        "sp|P1|A_YEAST again",
    ])
    summary = summarize_database(path)
    assert (summary.proteins, summary.decoys, summary.contaminants, summary.duplicates) == (5, 1, 1, 1)
    assert summary.organisms == {"Saccharomyces cerevisiae": 2, "Homo sapiens": 1}
    assert "5 proteins" in summary.text()
    ids = [n.id for n in check([], ORBITRAP_HCD, summary)]
    assert ids == ["database-duplicates"]


def test_a_large_or_empty_database_is_flagged(tmp_path, monkeypatch):
    from src.common import search_checks
    summary = summarize_database(_fasta(tmp_path / "db.fasta", [f"sp|P{i}|X" for i in range(50)]))
    monkeypatch.setattr(search_checks, "LARGE_DATABASE", 100)
    monkeypatch.setattr(search_checks, "LARGE_OPEN_DATABASE", 10)
    assert check([], ORBITRAP_HCD, summary) == []
    assert [n.id for n in check([], ORBITRAP_HCD, summary, open_search=True)] == ["database-large"]
    empty = summarize_database(_fasta(tmp_path / "empty.fasta", []))
    assert [n.id for n in check([], ORBITRAP_HCD, empty)] == ["database-empty"]
    missing = summarize_database(tmp_path / "missing.fasta")
    assert not missing.readable and [n.id for n in check([], ORBITRAP_HCD, missing)] == ["database-empty"]


def test_the_fingerprint_changes_with_what_the_notice_is_about():
    runs = [RunDetails("a", analyzer=instrument.ION_TRAP)]
    first = check(runs, ORBITRAP_HCD, None)[0].fingerprint
    second = check(runs, {**ORBITRAP_HCD, "Search:fragment:mass_tolerance": 10.0}, None)[0].fingerprint
    assert first != second
