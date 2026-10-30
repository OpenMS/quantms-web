"""Tests for the cRAP contaminants: database preparation and flagging in results."""
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pyopenms as poms

from src.common.contaminants import (
    CONTAMINANT_PREFIX,
    CRAP_FASTA,
    append_fasta,
    contaminant_mask,
    drop_contaminants,
    is_contaminant,
    write_contaminants,
)


def _fasta(path, entries):
    out = []
    for identifier, sequence in entries:
        e = poms.FASTAEntry()
        e.identifier = identifier
        e.sequence = sequence
        out.append(e)
    poms.FASTAFile().store(str(path), out)


def _ids(path):
    entries = []
    poms.FASTAFile().load(str(path), entries)
    return [e.identifier for e in entries]


def test_write_contaminants_prefixes_every_crap_entry(tmp_path):
    out = tmp_path / "con.fasta"
    n, skipped = write_contaminants(out)
    ids = _ids(out)
    assert n == len(ids) == len(_ids(CRAP_FASTA)) == 116
    assert skipped == 0
    assert all(i.startswith(CONTAMINANT_PREFIX) for i in ids)
    # An already prefixed entry is not prefixed twice
    write_contaminants(tmp_path / "again.fasta", crap_path=out)
    assert _ids(tmp_path / "again.fasta") == ids


def test_contaminants_already_in_the_database_stay_targets(tmp_path):
    crap = []
    poms.FASTAFile().load(str(CRAP_FASTA), crap)
    # One by sequence (under another name), one by UniProt entry name (another sequence)
    _fasta(tmp_path / "targets.fasta", [
        ("sp|Q0|MYPROT_BOVIN", crap[0].sequence),
        ("sp|P02768|ALBU_HUMAN", "PEPTIDEK"),
    ])
    assert any("ALBU_HUMAN" in i for i in _ids(CRAP_FASTA))
    out = tmp_path / "con.fasta"
    n, skipped = write_contaminants(out, targets=tmp_path / "targets.fasta")
    assert (n, skipped) == (114, 2)
    kept = _load(out)
    assert crap[0].sequence not in {e.sequence for e in kept}
    assert not any("ALBU_HUMAN" in e.identifier for e in kept)


def test_unchanged_contaminants_are_not_rewritten(tmp_path):
    out = tmp_path / "con.fasta"
    write_contaminants(out)
    before = out.stat().st_mtime_ns
    write_contaminants(out)
    assert out.stat().st_mtime_ns == before


def _load(path):
    entries = []
    poms.FASTAFile().load(str(path), entries)
    return entries


def test_append_fasta_keeps_order(tmp_path):
    _fasta(tmp_path / "a.fasta", [("sp|P1|A", "PEPTIDEK")])
    _fasta(tmp_path / "b.fasta", [("CON__sp|K1|B", "ELVISK"), ("rev_CON__sp|K1|B", "KSIVLE")])
    append_fasta([tmp_path / "a.fasta", tmp_path / "b.fasta"], tmp_path / "out" / "ab.fasta")
    assert _ids(tmp_path / "out" / "ab.fasta") == ["sp|P1|A", "CON__sp|K1|B", "rev_CON__sp|K1|B"]


def test_is_contaminant():
    assert is_contaminant("CON__sp|ALBU_BOVIN|")
    assert is_contaminant("sp|P1|A;CON__sp|K2C1_HUMAN|")
    assert is_contaminant("sp|P1|A,CON__X")
    assert is_contaminant("sp|P1|A/CON__X")
    assert not is_contaminant("sp|P1|A")
    assert not is_contaminant("rev_CON__X")  # a decoy, not a contaminant
    assert not is_contaminant(None)
    assert not is_contaminant(float("nan"))
    mask = contaminant_mask(pd.Series(["CON__X", "P1", None]))
    assert mask.tolist() == [True, False, False]
    df = pd.DataFrame({"ProteinName": ["CON__X", "P1"], "s1": [1.0, 2.0]})
    assert drop_contaminants(df, "ProteinName")["ProteinName"].tolist() == ["P1"]


def _workflow(tmp_path, params):
    """A WorkflowTest with only what prepare_search_database touches."""
    from src.WorkflowTest import WorkflowTest

    calls = []

    def run_topp(tool, io, custom_params={}, tool_instance_name=None):
        calls.append((tool, io, custom_params, tool_instance_name))
        # Stand in for DecoyDatabase: the targets, then each entry reversed with the decoy string
        decoy = custom_params.get("decoy_string", "rev_")
        inputs = io["in"][0] if isinstance(io["in"][0], list) else io["in"]
        entries = []
        for part in inputs:
            part_entries = []
            poms.FASTAFile().load(str(part), part_entries)
            entries.extend(part_entries)
        decoys = []
        for e in entries:
            d = poms.FASTAEntry()
            d.identifier = decoy + e.identifier
            d.sequence = e.sequence[::-1]
            decoys.append(d)
        poms.FASTAFile().store(io["out"][0], entries + decoys)
        return True

    wf = object.__new__(WorkflowTest)
    wf.params = params
    wf.workflow_dir = tmp_path / "topp-workflow"
    wf.results_dir = wf.workflow_dir / "results"
    wf.logger = SimpleNamespace(log=lambda *a, **k: None)
    wf.executor = SimpleNamespace(run_topp=run_topp)
    return wf, calls


def test_generated_decoys_cover_the_contaminants(tmp_path):
    fasta = tmp_path / "db.fasta"
    _fasta(fasta, [("sp|P1|A", "PEPTIDEK")])
    wf, calls = _workflow(tmp_path, {"generate-decoys": True})

    database, decoy_string = wf.prepare_search_database(fasta)

    assert decoy_string == "rev_"
    assert len(calls) == 1 and len(calls[0][1]["in"][0]) == 2  # FASTA and contaminants in one call
    ids = _ids(database)
    assert "sp|P1|A" in ids and "rev_sp|P1|A" in ids
    assert sum(i.startswith(CONTAMINANT_PREFIX) for i in ids) == 116
    assert sum(i.startswith("rev_" + CONTAMINANT_PREFIX) for i in ids) == 116
    # The user's FASTA is left alone
    assert _ids(fasta) == ["sp|P1|A"]


def test_fasta_with_decoys_gets_contaminants_in_its_decoy_style(tmp_path):
    fasta = tmp_path / "db.fasta"
    _fasta(fasta, [("sp|P1|A", "PEPTIDEK"), ("DECOY_sp|P1|A", "KEDITPEP")])
    params = {"generate-decoys": False, "CometAdapter": {"PeptideIndexing:decoy_string": "DECOY_"}}
    wf, calls = _workflow(tmp_path, params)

    database, decoy_string = wf.prepare_search_database(fasta)

    assert decoy_string == "DECOY_"
    assert calls[0][2]["decoy_string"] == "DECOY_"
    assert calls[0][3] == "DecoyDatabase-contaminants"
    ids = _ids(database)
    assert ids[:2] == ["sp|P1|A", "DECOY_sp|P1|A"]
    assert sum(i.startswith("DECOY_" + CONTAMINANT_PREFIX) for i in ids) == 116


def test_contaminants_can_be_switched_off(tmp_path):
    fasta = tmp_path / "db.fasta"
    _fasta(fasta, [("sp|P1|A", "PEPTIDEK"), ("rev_sp|P1|A", "KEDITPEP")])
    wf, calls = _workflow(tmp_path, {"generate-decoys": False, "add-contaminants": False})

    database, _ = wf.prepare_search_database(fasta)

    assert database == fasta and calls == []


def test_abundance_data_leaves_contaminants_out_by_default(tmp_path):
    from src.common.results_helpers import get_abundance_data

    workflow_dir = tmp_path / "topp-workflow"
    quant_dir = workflow_dir / "results" / "quant_results"
    quant_dir.mkdir(parents=True)
    (workflow_dir / "params.json").write_text(json.dumps({"analysis-mode": "LFQ"}))
    rows = []
    for protein in ("sp|P1|A", "CON__sp|ALBU_BOVIN|", "sp|P2|B;CON__sp|TRYP_PIG|"):
        for sample in ("s1", "s2"):
            rows.append({"ProteinName": protein, "PeptideSequence": "PEPTIDEK",
                         "Reference": f"{sample}.mzML", "Intensity": 1000.0})
    pd.DataFrame(rows).to_csv(quant_dir / "quant.csv", index=False)

    pivot_df, expr_df, _ = get_abundance_data(Path(tmp_path))
    assert pivot_df["ProteinName"].tolist() == ["sp|P1|A"]
    assert expr_df.index.tolist() == ["sp|P1|A"]

    pivot_all, expr_all, _ = get_abundance_data(Path(tmp_path), include_contaminants=True)
    assert len(pivot_all) == len(expr_all) == 3
