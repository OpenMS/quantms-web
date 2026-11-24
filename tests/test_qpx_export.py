"""Tests for src/common/qpx_export.py (QPX psm and provenance views)."""
# Sorted after test_parameter_defaults on purpose: writing PSMs imports
# src.common.results_helpers, which pulls in src.workflow.ParameterManager.
import json

import pyarrow.parquet as pq
import pytest
import pyopenms as poms

from src.common.qpx_export import (
    C13_DELTA,
    PROVENANCE_SCHEMA,
    PSM_SCHEMA,
    derive_id,
    logged_commands,
    mass_error_ppm,
    missed_cleavages,
    parse_command,
    provenance_rows,
    provenance_table,
    scan_number,
    write_psms,
    write_table,
    footer,
)

# Required fields of QPX 1.2's psm view (psm.yaml) and whether they may be null
PSM_REQUIRED = {
    "psm_id": False, "sequence": False, "peptidoform": False, "charge": False, "is_decoy": False,
    "calculated_mz": True, "observed_mz": False, "run_file_name": False, "scan": False,
}


def test_schemas_follow_qpx():
    for name, nullable in PSM_REQUIRED.items():
        assert PSM_SCHEMA.field(name).nullable is nullable, name
    assert str(PSM_SCHEMA.field("charge").type) == "int16"
    assert str(PSM_SCHEMA.field("scan").type) == "list<item: int32>"
    assert str(PSM_SCHEMA.field("posterior_error_probability").type) == "double"
    mod = PSM_SCHEMA.field("modifications").type.value_type
    assert [f.name for f in mod] == ["name", "accession", "positions"]
    assert [f.name for f in PROVENANCE_SCHEMA] == [
        "step_order", "step_category", "step_name", "tool_name", "tool_version",
        "tool_uri", "parameters", "config", "output_views",
    ]


def test_scan_numbers_and_missed_cleavages():
    assert scan_number("controllerType=0 controllerNumber=1 scan=52580") == 52580
    assert scan_number("index=12") == 12
    assert scan_number("17") == 17
    assert scan_number("spectrum=x") is None
    assert missed_cleavages("PEPKTIDER") == 1
    assert missed_cleavages("PEPKPTIDER") == 0  # K before P is not cleaved


def test_derive_id_is_stable_and_order_sensitive():
    a = derive_id(["PEPTIDEK", 2, "run1", [5]])
    assert a == derive_id(["PEPTIDEK", 2, "run1", [5]])
    assert a != derive_id(["PEPTIDEK", 2, "run1", [6]])
    assert -(2 ** 63) <= a < 2 ** 63


LOG = """STARTING WORKFLOW

Running command:
CometAdapter -in /ws/topp-workflow/input-files/mzML-files/a.mzML -out /ws/topp-workflow/runs/1/results/a.idXML -precursor_mass_tolerance -5 -threads 8 -threads 1
Waiting for command to finish...

Process finished:
CometAdapter -in /ws/topp-workflow/input-files/mzML-files/a.mzML -out /ws/topp-workflow/runs/1/results/a.idXML -precursor_mass_tolerance -5 -threads 8 -threads 1
Total time to run command: 1 seconds

Process finished:
CometAdapter -in /ws/topp-workflow/input-files/mzML-files/b.mzML -out /ws/topp-workflow/runs/1/results/b.idXML -precursor_mass_tolerance -5 -threads 8 -threads 1

Process finished:
IDFilter -in x.idXML -out y.idXML -score:psm 0.01
"""


def test_commands_are_read_from_the_log():
    commands = logged_commands(LOG)
    assert [c[0] for c in commands] == ["CometAdapter", "CometAdapter", "IDFilter"]
    known = {"in", "out", "precursor_mass_tolerance", "threads"}
    params = parse_command(commands[0], known)
    assert params["precursor_mass_tolerance"] == "-5"  # a negative value is not a parameter
    assert params["threads"] == "1"  # given twice: the tool takes the last


def fake_ini(tool):
    defaults = {
        "CometAdapter": {"in": "", "out": "", "precursor_mass_tolerance": 10.0, "threads": 1, "enzyme": "Trypsin",
                         "database": ""},
        "IDFilter": {"in": "", "out": "", "score:psm": 0.0},
    }[tool]
    return "3.5.0", defaults


def test_provenance_lists_every_tool_and_parameter():
    params = {"analysis-mode": "LFQ", "fasta-file": "/ws/topp-workflow/input-files/fasta-file/db.fasta",
              "psm-fdr-percent": 1.0, "mzML-files": ["/x/a.mzML", "/x/b.mzML"]}
    rows = provenance_rows(LOG, "/ws/topp-workflow", params, "1.2.1", ini=fake_ini, engine=lambda t, p: "2026.01" if t == "CometAdapter" else None)
    assert [r["tool_name"] for r in rows] == ["OpenDDA", "OpenMS/CometAdapter", "OpenMS/IDFilter", "OpenDDA"]
    assert [r["step_order"] for r in rows] == [1, 2, 3, 4]
    comet = {p["key"]: p["value"] for p in rows[1]["parameters"]}
    # one row per tool; per-file values list every file, paths relative to the workflow folder
    assert comet["in"] == "input-files/mzML-files/a.mzML; input-files/mzML-files/b.mzML"
    assert comet["engine_version"] == "2026.01"
    config = json.loads(rows[1]["config"])
    assert config["enzyme"] == "Trypsin"  # every parameter, defaults included
    assert config["precursor_mass_tolerance"] == "-5"
    workflow = {p["key"]: p["value"] for p in rows[0]["parameters"]}
    assert workflow["fasta_file"] == "input-files/fasta-file/db.fasta"
    assert workflow["mzML_files"] == "a.mzML b.mzML"

    table = provenance_table(rows)
    comet_rows = table[table["tool"] == "OpenMS/CometAdapter"]
    assert set(comet_rows["parameter"]) >= {"enzyme", "precursor_mass_tolerance", "engine_version"}
    assert comet_rows.set_index("parameter").loc["enzyme", "set by"] == "tool default"
    assert comet_rows.set_index("parameter").loc["threads", "set by"] == "workflow"


def _idxml(path, isotope_errors=None):
    prot = poms.ProteinIdentification()
    prot.setIdentifier("run")
    hits = []
    for acc in ("P1", "rev_P1"):
        protein_hit = poms.ProteinHit()
        protein_hit.setAccession(acc)
        hits.append(protein_hit)
    prot.setHits(hits)
    peptides = poms.PeptideIdentificationList() if hasattr(poms, "PeptideIdentificationList") else []
    for i, (seq, td, acc) in enumerate([("PEPTM(Oxidation)IDEK", "target", "P1"), ("KEDITPEP", "decoy", "rev_P1")]):
        isotope_error = isotope_errors[i] if isotope_errors else 0
        pep = poms.PeptideIdentification()
        pep.setIdentifier("run")
        pep.setScoreType("Posterior Error Probability")
        pep.setHigherScoreBetter(False)
        # a hit matched k isotopes off (k = -isotope_error) is observed k * 13C/z above its m/z
        pep.setMZ(poms.AASequence.fromString(seq).getMZ(2) - isotope_error * C13_DELTA / 2 if isotope_errors else 500.0 + i)
        pep.setRT(100.0 * (i + 1))
        pep.setMetaValue("spectrum_reference", f"controllerType=0 controllerNumber=1 scan={10 + i}")
        hit = poms.PeptideHit()
        hit.setSequence(poms.AASequence.fromString(seq))
        hit.setCharge(2)
        hit.setScore(0.01)
        hit.setMetaValue("target_decoy", td)
        hit.setMetaValue("MS:1001491", 0.002)
        if isotope_errors:
            hit.setMetaValue("isotope_error", isotope_error)
        ev = poms.PeptideEvidence()
        ev.setProteinAccession(acc)
        hit.setPeptideEvidences([ev])
        pep.setHits([hit])
        peptides.push_back(pep) if hasattr(peptides, "push_back") else peptides.append(pep)
    poms.IdXMLFile().store(str(path), [prot], peptides)


def test_psm_view_from_an_idxml(tmp_path):
    idxml = tmp_path / "sample1_filter.idXML"
    _idxml(idxml)
    out = tmp_path / "qpx" / "run.psm.parquet"
    n = write_psms([idxml], ["sample1.mzML"], out, to_proforma=lambda s: s.replace("(Oxidation)", "[Oxidation]"))
    assert n == 2
    table = pq.read_table(out)
    assert table.schema.equals(PSM_SCHEMA, check_metadata=False)
    meta = table.schema.metadata
    assert meta[b"qpx_version"] == b"1.2" and meta[b"file_type"] == b"psm_file" and meta[b"scan_format"] == b"scan"
    rows = table.to_pylist()
    target, decoy = rows
    assert target["run_file_name"] == "sample1"
    assert target["scan"] == [10]
    assert target["peptidoform"] == "PEPTM[Oxidation]IDEK"
    assert target["sequence"] == "PEPTMIDEK"
    assert target["modifications"][0]["accession"] == "UNIMOD:35"
    assert target["modifications"][0]["positions"][0] == {"position": 5, "amino_acid": "M", "scores": None}
    assert abs(target["posterior_error_probability"] - 0.01) < 1e-6
    assert target["protein_accessions"] == ["P1"]
    assert [s["score_name"] for s in target["additional_scores"]] == ["percolator q-value"]
    assert abs(target["calculated_mz"] - poms.AASequence.fromString("PEPTM(Oxidation)IDEK").getMZ(2)) < 1e-3
    assert not target["is_decoy"] and decoy["is_decoy"]
    assert target["psm_id"] == derive_id([target["peptidoform"], 2, "sample1", [10]])


def test_mass_error_is_to_the_matched_isotope(tmp_path):
    assert mass_error_ppm(500.005, 500.0, 2) == pytest.approx(10.0)
    # one isotope above at charge 2, reported by the engine as isotope_error -1
    assert mass_error_ppm(500.0 + C13_DELTA / 2 + 0.005, 500.0, 2, -1) == pytest.approx(10.0)
    assert mass_error_ppm(500.0, None, 2) is None

    idxml = tmp_path / "sample1_filter.idXML"
    _idxml(idxml, isotope_errors=[-1, 0])
    write_psms([idxml], ["sample1.mzML"], tmp_path / "run.psm.parquet", to_proforma=lambda s: s)
    errors = pq.read_table(tmp_path / "run.psm.parquet").column("mass_error_ppm").to_pylist()
    assert errors == [pytest.approx(0.0, abs=0.01), pytest.approx(0.0, abs=0.01)]


def test_provenance_view_round_trip(tmp_path):
    rows = provenance_rows(LOG, "/ws/topp-workflow", {}, "1.2.1", ini=fake_ini, engine=lambda t, p: None)
    path = write_table(rows, PROVENANCE_SCHEMA, tmp_path / "run.provenance.parquet", footer("provenance_file"))
    back = pq.read_table(path)
    assert back.schema.equals(PROVENANCE_SCHEMA, check_metadata=False)
    assert back.num_rows == len(rows)


# ----------------------------------------------------------------------------- pg and run views

def _msstats():
    import pandas as pd

    return pd.DataFrame([
        # ProteinName, PeptideSequence, Charge, ..., Reference, Intensity
        ("P1;P2", "PEPTIDEK", 2, "A", 1, 1, "a_F1.mzML", 100.0),
        ("P1;P2", "PEPM(Oxidation)TIDEK", 2, "A", 1, 1, "a_F1.mzML", 50.0),
        ("P1;P2", "PEPMTIDEK", 3, "A", 1, 1, "a_F2.mzML", 25.0),
        ("P1;P2", "PEPTIDEK", 2, "B", 2, 2, "b.mzML", 10.0),
        ("DECOY_P3", "KEDITPEP", 2, "B", 2, 2, "b.mzML", 5.0),
        ("CON__P4", "SAMPLER", 2, "B", 2, 2, "b.mzML", 7.0),
    ], columns=["ProteinName", "PeptideSequence", "PrecursorCharge", "Condition", "BioReplicate", "Run",
                "Reference", "Intensity"])


FRACTIONS = {
    "a_F1.mzML": {"sample": "A_1", "condition": "A"},
    "a_F2.mzML": {"sample": "A_1", "condition": "A"},
    "b.mzML": {"sample": "B_1", "condition": "B"},
}


def pa_table(rows, schema):
    import pyarrow as pa

    return pa.Table.from_pylist(rows, schema=schema)


def test_pg_rows_sum_a_group_per_sample_and_group_fractions():
    from src.common.qpx_export import PG_SCHEMA, pg_rows

    rows = pg_rows(_msstats(), FRACTIONS)
    by = {(tuple(r["pg_accessions"]), tuple(r["grouped_runs"])): r for r in rows}
    a = by[(("P1", "P2"), ("a_F1", "a_F2"))]
    assert a["intensity"] == pytest.approx(175.0)
    assert a["anchor_protein"] == "P1" and a["label"] == "LFQ"
    # PEPTIDEK, PEPMTIDEK and its oxidized form: three peptides, two sequences
    assert a["peptide_counts"] == {"unique_sequences": 2, "total_sequences": 3}
    assert a["pg_id"] == derive_id([["P1", "P2"], ["a_F1", "a_F2"], "LFQ"])
    assert a["pg_group_id"] == by[(("P1", "P2"), ("b",))]["pg_group_id"]
    assert a["pg_id"] != by[(("P1", "P2"), ("b",))]["pg_id"]
    assert by[(("DECOY_P3",), ("b",))]["is_decoy"] is True
    assert by[(("CON__P4",), ("b",))]["contaminant"] is True
    assert not a["is_decoy"] and not a["contaminant"]
    assert pa_table(rows, PG_SCHEMA).num_rows == len(rows)


def test_pg_rows_without_fractions_use_one_file_per_sample():
    from src.common.qpx_export import pg_rows

    rows = pg_rows(_msstats(), None)
    assert {tuple(r["grouped_runs"]) for r in rows if r["anchor_protein"] == "P1"} == {("a_F1",), ("a_F2",), ("b",)}


def test_run_rows_carry_sample_replicate_fraction_and_enzyme():
    from src.common.qpx_export import run_rows

    params = {
        "mzML-files": ["/ws/in/a_F1.mzML", "/ws/in/a_F2.mzML", "/ws/in/b.mzML"],
        "_defaults": {"ProSE:1:Search:enzyme": "Trypsin/P"},
    }
    rows = {r["run_file_name"]: r for r in run_rows(params, _msstats(), FRACTIONS)}
    assert set(rows) == {"a_F1", "a_F2", "b"}
    assert rows["a_F2"]["samples"][0]["sample_accession"] == "A_1"
    assert rows["a_F2"]["fraction"] == "2" and rows["b"]["fraction"] == "1"
    assert rows["b"]["samples"][0]["biological_replicate"] == 2
    assert rows["a_F1"]["file_name"] == "a_F1.mzML"
    assert rows["a_F1"]["enzymes"] == ["Trypsin/P"] and rows["a_F1"]["acquisition_method"] == "DDA"


def test_export_views_write_qpx_pg_and_run_files(tmp_path):
    from src.common.qpx_export import PG_SCHEMA, RUN_SCHEMA, export_views

    quant = tmp_path / "quant_results"
    quant.mkdir()
    _msstats().to_csv(quant / "openms_msstats.csv", index=False)
    out = export_views(tmp_path, {"mzML-files": ["a_F1.mzML", "a_F2.mzML", "b.mzML"]})
    assert out["protein_groups"] == 5  # P1;P2 in each of three files, the decoy, the contaminant
    pg = pq.read_table(out["pg"])
    assert pg.schema.equals(PG_SCHEMA, check_metadata=False)
    meta = {k.decode(): v.decode() for k, v in pg.schema.metadata.items()}
    assert meta["file_type"] == "pg_file" and meta["identity_composite"] == "pg_accessions,grouped_runs,label"
    run = pq.read_table(out["run"])
    assert run.schema.equals(RUN_SCHEMA, check_metadata=False) and run.num_rows == 3
    # An isobaric run without its protein table gets the run view only
    iso = export_views(tmp_path / "iso", {"analysis-mode": "TMT", "mzML-files": ["x.mzML"]})
    assert "pg" not in iso and pq.read_table(iso["run"]).num_rows == 1


def test_run_enzyme_comes_from_the_search_step_of_the_provenance():
    from src.common.qpx_export import run_rows

    steps = [{"step_name": "sequence_search", "parameters": [{"key": "database", "value": "db.fasta"}],
              "config": json.dumps({"Search:enzyme": "Lys-C"})}]
    rows = run_rows({"mzML-files": ["a.mzML"]}, None, None, steps)
    assert rows[0]["enzymes"] == ["Lys-C"] and rows[0]["samples"][0]["sample_accession"] == "a"


ISOBARIC_TABLE = (
    "# Protein abundances computed from file 'merged.consensusXML'\n"
    "protein\tn_proteins\tprotein_score\tn_peptides\tabundance_sample1[unknown]\t"
    "abundance_sample2[unknown]\tabundance_sample3[unknown]\tratio_sample2[unknown]/sample1[unknown]\n"
    "P1/P2\t2\t0.99\t5\t100\t200\t0\t2\n"
    "DECOY_P3\t1\t0.5\t1\t10\t20\t30\t2\n"
    "CON__P4\t1\t0.9\t2\t1\tnan\t3\t1\n"
)
TMT_PARAMS = {
    "analysis-mode": "TMT",
    "mzML-files": ["set_F1.mzML", "set_F2.mzML"],
    "IsobaricAnalyzer-TMT": {"type": "tmt11plex"},
}


def _isobaric_run(tmp_path):
    quant = tmp_path / "quant_results"
    quant.mkdir(parents=True)
    (quant / "openms_design_protein_openms.csv").write_text(ISOBARIC_TABLE)
    return tmp_path


def test_isobaric_pg_rows_one_per_group_and_channel(tmp_path):
    from src.common.qpx_export import derive_id, isobaric_pg_rows, read_isobaric_proteins

    proteins = read_isobaric_proteins(_isobaric_run(tmp_path))
    rows = isobaric_pg_rows(proteins, "tmt11plex", ["set_F2", "set_F1"])
    assert len(rows) == 9  # three groups, three channels; the ratio column is not a channel
    first = rows[:3]
    assert [r["label"] for r in first] == ["TMT126", "TMT127N", "TMT127C"]
    assert [r["intensity"] for r in first] == [100.0, 200.0, 0.0]
    assert first[0]["pg_accessions"] == ["P1", "P2"] and first[0]["anchor_protein"] == "P1"
    assert first[0]["grouped_runs"] == ["set_F1", "set_F2"]
    assert first[0]["pg_id"] == derive_id([["P1", "P2"], ["set_F1", "set_F2"], "TMT126"])
    assert len({r["pg_id"] for r in rows}) == 9 and len({r["pg_group_id"] for r in first}) == 1
    assert first[0]["peptide_counts"] == {"unique_sequences": 5, "total_sequences": 5}
    assert rows[3]["is_decoy"] and not rows[3]["contaminant"]
    assert rows[6]["contaminant"] and rows[7]["intensity"] is None


def test_itraq_channels_are_labelled_itraq():
    import pandas as pd
    from src.common.qpx_export import isobaric_channels

    cols = ["protein", "abundance_sample2[unknown]", "abundance_sample1[unknown]"]
    assert isobaric_channels(pd.Index(cols), "itraq4plex") == [
        ("abundance_sample1[unknown]", 1, "iTRAQ114"),
        ("abundance_sample2[unknown]", 2, "iTRAQ115"),
    ]


def test_export_views_write_an_isobaric_pg_and_channel_samples(tmp_path):
    from src.common.qpx_export import PG_SCHEMA, export_views, views_missing

    results = _isobaric_run(tmp_path)
    assert views_missing(results, TMT_PARAMS)
    out = export_views(results, TMT_PARAMS)
    assert out["protein_groups"] == 9
    pg = pq.read_table(out["pg"])
    assert pg.schema.equals(PG_SCHEMA, check_metadata=False)
    runs = pq.read_table(out["run"]).to_pylist()
    assert [r["run_file_name"] for r in runs] == ["set_F1", "set_F2"]
    assert [(s["sample_accession"], s["label"]) for s in runs[0]["samples"]] == [
        ("sample1", "TMT126"), ("sample2", "TMT127N"), ("sample3", "TMT127C"),
    ]
    assert not views_missing(results, TMT_PARAMS)


def test_isobaric_run_view_without_pg_is_completed(tmp_path):
    from src.common.qpx_export import RUN_SCHEMA, SAMPLE_SCHEMA, view_paths, views_missing, write_table

    results = _isobaric_run(tmp_path)
    # A run exported before the isobaric pg view: the run and sample views are there, the pg view not
    write_table([], RUN_SCHEMA, view_paths(results)["run"], {})
    write_table([], SAMPLE_SCHEMA, view_paths(results)["sample"], {})
    assert views_missing(results, TMT_PARAMS)
    # Without the protein table there is nothing to add
    bare = tmp_path / "bare"
    write_table([], RUN_SCHEMA, view_paths(bare)["run"], {})
    write_table([], SAMPLE_SCHEMA, view_paths(bare)["sample"], {})
    assert not views_missing(bare, TMT_PARAMS)


# ----------------------------------------------------------------------------- sample view

FASTA = (
    ">sp|P1|A_HUMAN Protein A OS=Homo sapiens OX=9606 GN=A PE=1 SV=1\nPEPTIDEK\n"
    ">sp|P2|B_HUMAN Protein B OS=Homo sapiens OX=9606 GN=B PE=1 SV=1\nPEPTIDER\n"
    ">sp|P3|C_YEAST Protein C OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 PE=1 SV=1\nKEDITPEP\n"
    ">DECOY_sp|P1|A_HUMAN\nKEDITPEP\n"
)
SDRF = (
    "source name\tcharacteristics[organism]\tcharacteristics[organism part]\tcharacteristics[disease]\t"
    "characteristics[sex]\tcomment[label]\tcomment[data file]\n"
    "patient 1\tHomo sapiens\tliver\thepatocellular carcinoma\tfemale\tlabel free sample\ta_F1.raw\n"
    "patient 1\tHomo sapiens\tliver\thepatocellular carcinoma\tfemale\tlabel free sample\ta_F2.raw\n"
    "patient 2\tHomo sapiens\tliver\tnormal\tnot available\tlabel free sample\tb.raw\n"
)


def test_the_database_names_the_organism(tmp_path):
    from src.common.qpx_export import fasta_organism

    fasta = tmp_path / "db.fasta"
    fasta.write_text(FASTA)
    assert fasta_organism(fasta) == "Homo sapiens"
    (tmp_path / "plain.fasta").write_text(">P1\nPEPTIDEK\n")
    assert fasta_organism(tmp_path / "plain.fasta") is None
    assert fasta_organism(tmp_path / "missing.fasta") is None


def test_sample_rows_take_condition_replicate_and_sdrf_characteristics(tmp_path):
    from src.common.qpx_export import SAMPLE_SCHEMA, run_rows, sample_rows, sdrf_samples

    sdrf = tmp_path / "PXD1.sdrf.tsv"
    sdrf.write_text(SDRF)
    params = {"mzML-files": ["a_F1.mzML", "a_F2.mzML", "b.mzML"],
              "mzML-group-a_F1.mzML": "A", "mzML-group-b.mzML": "B"}
    runs = run_rows(params, _msstats(), FRACTIONS)
    rows = {r["sample_accession"]: r for r in sample_rows(runs, params, FRACTIONS, sdrf_samples([sdrf]), "Mus musculus")}
    assert set(rows) == {"A_1", "B_1"}  # the two fractions of A_1 are one sample
    a, b = rows["A_1"], rows["B_1"]
    assert a["organism"] == "Homo sapiens" and a["organism_part"] == "liver"
    assert a["disease"] == "hepatocellular carcinoma" and a["sex"] == "female"
    assert a["condition"] == "A" and a["biological_replicate"] == 1 and a["source_name"] == "patient 1"
    assert b["condition"] == "B" and b["biological_replicate"] == 2
    assert b["sex"] is None  # "not available" in the SDRF
    assert pa_table(list(rows.values()), SAMPLE_SCHEMA).num_rows == 2


def test_sample_rows_without_sdrf_take_the_database_organism():
    from src.common.qpx_export import NOT_AVAILABLE, run_rows, sample_rows

    params = {"mzML-files": ["a.mzML", "b.mzML"], "mzML-group-a.mzML": "ctrl"}
    rows = sample_rows(run_rows(params), params, organism="Homo sapiens")
    assert [r["sample_accession"] for r in rows] == ["a", "b"]
    assert rows[0]["organism"] == "Homo sapiens" and rows[0]["organism_part"] == NOT_AVAILABLE
    assert rows[0]["condition"] == "ctrl" and rows[1]["condition"] is None
    assert sample_rows(run_rows(params), params)[0]["organism"] == NOT_AVAILABLE


def test_isobaric_samples_are_the_channels_with_their_groups(tmp_path):
    from src.common.qpx_export import SAMPLE_SCHEMA, export_views

    results = _isobaric_run(tmp_path / "wf" / "runs" / "r1" / "results")
    sdrf_dir = tmp_path / "wf" / "pride-import" / "sdrf"
    sdrf_dir.mkdir(parents=True)
    (sdrf_dir / "PXD2.sdrf.tsv").write_text(
        "source name\tcharacteristics[organism]\tcharacteristics[organism part]\tcomment[label]\tcomment[data file]\n"
        "pool\tMus musculus\tbrain\tTMT126\tset_F1.raw\n"
        "treated\tMus musculus\tbrain\tTMT127N\tset_F1.raw\n"
    )
    params = {**TMT_PARAMS, "TMT-group-sample1": "ref", "TMT-group-sample2": "drug"}
    out = export_views(results, params)
    table = pq.read_table(out["sample"])
    assert table.schema.equals(SAMPLE_SCHEMA, check_metadata=False)
    meta = {k.decode(): v.decode() for k, v in table.schema.metadata.items()}
    assert meta["file_type"] == "sample_file" and meta["primary_key"] == "sample_accession"
    rows = table.to_pylist()
    assert [r["sample_accession"] for r in rows] == ["sample1", "sample2", "sample3"]
    assert [r["condition"] for r in rows] == ["ref", "drug", None]
    assert [r["source_name"] for r in rows] == ["pool", "treated", None]
    assert rows[1]["organism"] == "Mus musculus" and rows[2]["organism"] == "not available"
