"""Verdicts on single PSMs, kept in results/curation.tsv."""

import pandas as pd
import pytest

from src.common.psm_curation import CURATION_FILE, load_verdicts, set_verdict, verdict_of, verdicts_tsv

PSM = {"file": "a.mzML", "scan": 1234, "sequence": "PEPTM(Oxidation)IDEK", "charge": 2, "mz": 547.27, "rt": 1802.5}


def test_no_file_means_no_verdicts(tmp_path):
    frame = load_verdicts(tmp_path)
    assert frame.empty
    assert verdict_of(frame, PSM) == (None, "")


def test_a_verdict_and_note_are_kept_and_read_back(tmp_path):
    set_verdict(tmp_path, PSM, "rejected", "few fragments matched")
    frame = load_verdicts(tmp_path)
    assert len(frame) == 1
    assert verdict_of(frame, PSM) == ("rejected", "few fragments matched")
    assert (tmp_path / CURATION_FILE).is_file()


def test_a_new_verdict_replaces_the_old_one(tmp_path):
    set_verdict(tmp_path, PSM, "rejected", "first look")
    set_verdict(tmp_path, PSM, "accepted", "b and y series complete")
    frame = load_verdicts(tmp_path)
    assert len(frame) == 1
    assert verdict_of(frame, PSM) == ("accepted", "b and y series complete")


def test_psms_differ_by_file_scan_peptide_and_charge(tmp_path):
    set_verdict(tmp_path, PSM, "accepted")
    for change in ({"file": "b.mzML"}, {"scan": 1235}, {"sequence": "PEPTMIDEK"}, {"charge": 3}):
        assert verdict_of(load_verdicts(tmp_path), {**PSM, **change}) == (None, "")
    set_verdict(tmp_path, {**PSM, "charge": 3}, "rejected")
    assert len(load_verdicts(tmp_path)) == 2


def test_clearing_verdict_and_note_takes_the_psm_off(tmp_path):
    set_verdict(tmp_path, PSM, "accepted", "fine")
    set_verdict(tmp_path, PSM, None, "")
    assert load_verdicts(tmp_path).empty


def test_a_note_alone_is_kept(tmp_path):
    set_verdict(tmp_path, PSM, None, "check the isotope pattern")
    assert verdict_of(load_verdicts(tmp_path), PSM) == (None, "check the isotope pattern")


def test_tabs_and_newlines_in_a_note_do_not_break_the_table(tmp_path):
    set_verdict(tmp_path, PSM, "rejected", "line one\nline\ttwo")
    assert verdict_of(load_verdicts(tmp_path), PSM) == ("rejected", "line one line two")


def test_an_unknown_verdict_is_refused(tmp_path):
    with pytest.raises(ValueError):
        set_verdict(tmp_path, PSM, "maybe")


def test_an_unreadable_file_reads_as_no_verdicts(tmp_path):
    (tmp_path / CURATION_FILE).write_bytes(b"\x00\x01 not a table \"")
    assert load_verdicts(tmp_path).empty


def test_the_download_lists_accepted_first(tmp_path):
    set_verdict(tmp_path, PSM, "rejected", "no")
    set_verdict(tmp_path, {**PSM, "scan": 99}, "accepted", "yes")
    lines = verdicts_tsv(load_verdicts(tmp_path)).decode().splitlines()
    assert lines[0].split("\t")[:4] == ["file", "scan", "sequence", "charge"]
    assert "accepted" in lines[1] and "rejected" in lines[2]
    table = pd.read_csv(tmp_path / CURATION_FILE, sep="\t")
    assert set(table.columns) >= {"verdict", "note", "updated", "mz", "rt"}
