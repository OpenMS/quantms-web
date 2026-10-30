import numpy as np
import pandas as pd

from src.common.protein_table import (
    build_tables,
    coverage_html,
    coverage_mask,
    gene_name,
    protein_name,
    read_fasta,
)


def test_coverage_mask_counts_every_occurrence():
    mask = coverage_mask("PEPTIDEKAAAPEPTIDEK", ["PEPTIDEK", "AA"])
    assert mask.tolist() == [True] * 8 + [True, True, True] + [True] * 8
    assert coverage_mask("ACDEFG", ["XYZ"]).sum() == 0


def test_uniprot_description_fields():
    desc = "Serum albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2"
    assert gene_name(desc) == "ALB"
    assert protein_name(desc) == "Serum albumin"
    assert gene_name("no gene here") == ""


def test_read_fasta_skips_decoys(tmp_path):
    fasta = tmp_path / "db.fasta"
    fasta.write_text(">sp|P1|A_HUMAN Protein A GN=A\nPEPTIDEKAAA\n>rev_sp|P1|A_HUMAN\nAAAKEDITPEP\n")
    entries = read_fasta(fasta)
    assert list(entries) == ["sp|P1|A_HUMAN"]
    assert entries["sp|P1|A_HUMAN"][1] == "PEPTIDEKAAA"


def test_build_tables():
    fasta = {
        "P1": ("Protein one GN=ONE", "PEPTIDEKSAMPLERAAAA"),
        "P2": ("Protein two GN=TWO", "SAMPLERGGGG"),
    }
    psms = pd.DataFrame({
        "file": ["a.mzML", "a.mzML", "b.mzML"],
        "scan": [1, 2, 3],
        "rt": [10.0, 20.0, 30.0],
        "mz": [400.0, 500.0, 400.0],
        "charge": [2, 2, 3],
        "peptide": ["PEPTIDEK", "SAMPLER", "PEPTIDEK"],
        "sequence": ["PEPTIDEK", "SAMPLER", "PEPTIDEK"],
        "score": [0.01, 0.02, 0.03],
        "q_value": [0.001, 0.002, 0.003],
        "proteins": [["P1"], ["P1", "P2"], ["P1"]],
    })
    proteins, peptides, psms_long, coverage = build_tables(psms, fasta, ["a.mzML", "b.mzML"])

    p1 = proteins.set_index("accession").loc["P1"]
    assert p1["peptides"] == 2
    assert p1["unique peptides"] == 1  # SAMPLER is shared with P2
    assert p1["PSMs"] == 3
    assert p1["found in"] == "2 of 2"
    assert p1["coverage (%)"] == round(100 * 15 / 19, 1)
    assert p1["PSMs a"] == 2 and p1["PSMs b"] == 1
    assert proteins.iloc[0]["accession"] == "P1"  # most PSMs first

    assert set(peptides.loc[peptides["accession"] == "P2", "peptide"]) == {"SAMPLER"}
    assert len(psms_long) == 4  # the shared PSM appears under both proteins
    assert coverage["P2"].tolist() == [True] * 7 + [False] * 4


def test_coverage_html_marks_covered_residues():
    out = coverage_html("ACDEFGHIKLM", np.array([True, True] + [False] * 9))
    assert out.count('<span class="cov">') == 2
    assert "ACDEFGHIKL M" in out.replace('<span class="cov">A</span>', "A").replace('<span class="cov">C</span>', "C")


def test_psm_ids_match_the_spectrum_viewer_numbering(tmp_path):
    """collect_psms' psm_id is the top hit's id_idx in parse_idxml (the viewer cache)."""
    import pyopenms as poms

    from src.common.protein_table import collect_psms
    from src.common.results_helpers import parse_idxml

    prot = poms.ProteinIdentification()
    prot.setIdentifier("run")
    protein_hits = []
    for acc in ("P1", "P2", "rev_P1"):
        protein_hit = poms.ProteinHit()
        protein_hit.setAccession(acc)
        protein_hits.append(protein_hit)
    prot.setHits(protein_hits)
    peptides = []
    # Two hits, then a decoy-only spectrum, then one hit: the ids skip every hit
    for scan, hits in [(5, [("PEPTIDEK", "P1"), ("PEPTLDEK", "P2")]), (6, [("KEDITPEP", "rev_P1")]), (9, [("SAMPLER", "P2")])]:
        pep = poms.PeptideIdentification()
        pep.setIdentifier("run")
        pep.setMetaValue("spectrum_reference", f"controllerType=0 controllerNumber=1 scan={scan}")
        out = []
        for seq, acc in hits:
            hit = poms.PeptideHit()
            hit.setSequence(poms.AASequence.fromString(seq))
            hit.setCharge(2)
            ev = poms.PeptideEvidence()
            ev.setProteinAccession(acc)
            hit.setPeptideEvidences([ev])
            out.append(hit)
        pep.setHits(out)
        peptides.append(pep)
    if hasattr(poms, "PeptideIdentificationList"):
        container = poms.PeptideIdentificationList()
        for pep in peptides:
            container.push_back(pep)
        peptides = container
    path = tmp_path / "runA_filter.idXML"
    poms.IdXMLFile().store(str(path), [prot], peptides)

    psms = collect_psms([path], ["runA.mzML"])
    viewer = parse_idxml(path, {"runA.mzML": 0})[0].to_pandas().set_index("id_idx")
    assert psms["psm_id"].tolist() == [0, 3]
    assert psms["run"].tolist() == ["runA_filter", "runA_filter"]
    for _, row in psms.iterrows():
        assert viewer.loc[row["psm_id"], "sequence"] == row["peptide"]
        assert viewer.loc[row["psm_id"], "scan_id"] == row["scan"]
