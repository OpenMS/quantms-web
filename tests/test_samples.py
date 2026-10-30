"""Tests for the per-sample counts (Samples page)."""
import pandas as pd
import pyopenms as poms

# src.common.samples is imported inside the tests, not at collection: importing
# it pulls in src.workflow.ParameterManager with the real streamlit, and
# test_parameter_defaults needs to import it with a mock.


def _idxml(path, hits):
    """hits: (sequence, protein accession, "target" or "decoy")."""
    prot = poms.ProteinIdentification()
    prot.setIdentifier("run")
    peptides = []
    for seq, accession, td in hits:
        pep = poms.PeptideIdentification()
        pep.setIdentifier("run")
        hit = poms.PeptideHit()
        hit.setSequence(poms.AASequence.fromString(seq))
        hit.setCharge(2)
        hit.setMetaValue("target_decoy", td)
        evidence = poms.PeptideEvidence()
        evidence.setProteinAccession(accession)
        hit.setPeptideEvidences([evidence])
        pep.setHits([hit])
        peptides.append(pep)
    protein_hits = []
    for accession in sorted({h[1] for h in hits}):
        protein_hit = poms.ProteinHit()
        protein_hit.setAccession(accession)
        protein_hits.append(protein_hit)
    prot.setHits(protein_hits)
    if hasattr(poms, "PeptideIdentificationList"):
        container = poms.PeptideIdentificationList()
        for pep in peptides:
            container.push_back(pep)
        peptides = container
    path.parent.mkdir(parents=True, exist_ok=True)
    poms.IdXMLFile().store(str(path), [prot], peptides)


def _quant_csv(path, rows):
    """rows: (protein, peptide, reference, intensity)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["ProteinName", "PeptideSequence", "Reference", "Intensity"]).to_csv(path, index=False)


def test_identification_counts_skip_decoys(tmp_path):
    from src.common.samples import identification_counts

    _idxml(tmp_path / "a.idXML", [
        ("PEPTIDE", "P1", "target"),
        ("PEPTIDE", "P1", "target"),
        ("PEPTIDEK", "P2", "target"),
        ("EDITPEP", "DECOY_P1", "decoy"),
    ])
    assert identification_counts(tmp_path / "a.idXML") == {"Target PSMs": 3, "Peptides": 2, "Proteins": 2}


def test_lfq_samples_group_by_condition(tmp_path):
    from src.common.samples import lfq_samples

    results = tmp_path / "results"
    _idxml(results / "psm_filter" / "ctrl1_filter.idXML", [("PEPTIDE", "P1", "target")])
    _idxml(results / "psm_filter" / "case1_filter.idXML", [
        ("PEPTIDE", "P1", "target"), ("PEPTIDEK", "P2", "target"),
    ])
    _idxml(results / "psm_filter" / "case2_filter.idXML", [("PEPTIDE", "P1", "target")])
    _quant_csv(results / "quant_results" / "quant.csv", [
        ("P1", "PEPTIDE", "ctrl1.mzML", 10.0),
        ("P1", "PEPTIDE", "case1.mzML", 5.0),
        ("P2", "PEPTIDEK", "case1.mzML", 3.0),
        ("P2", "PEPTIDEK", "case2.mzML", 0.0),
    ])
    files = ["in/ctrl1.mzML", "in/case1.mzML", "in/case2.mzML", "in/blank.mzML"]
    groups = {"ctrl1.mzML": "control", "case1.mzML": "case", "case2.mzML": "case"}

    df = lfq_samples(tmp_path, files, groups)

    # Conditions sorted, unassigned last; file order kept within a condition
    assert df["Sample"].tolist() == ["case1", "case2", "ctrl1", "blank"]
    assert df["Condition"].tolist() == ["case", "case", "control", "Unassigned"]
    case1 = df.iloc[0]
    assert (case1["Target PSMs"], case1["Peptides"], case1["Proteins"]) == (2, 2, 2)
    assert (case1["Quantified proteins"], case1["Quantified peptides"]) == (2, 2)
    # An intensity of 0 is not quantified
    assert pd.isna(df.iloc[1]["Quantified proteins"])
    # A file the run did not reach has empty cells
    assert df.iloc[3][["Target PSMs", "Quantified proteins"]].isna().all()


def test_channel_samples_and_condition_summary():
    from src.common.samples import channel_samples, condition_summary

    pivot = pd.DataFrame({
        "protein": ["P1", "P2", "P3"],
        "sample1[126]": [1.0, 2.0, 0.0],
        "sample2[127N]": [1.0, None, None],
        "sample3[127C]": [1.0, 1.0, 1.0],
    })
    groups = {"sample1[126]": "A", "sample2[127N]": "A", "sample3[127C]": "B"}
    df = channel_samples(pivot, ["sample1[126]", "sample2[127N]", "sample3[127C]"], groups)
    assert df["Quantified proteins"].tolist() == [2, 1, 3]
    assert df["Target PSMs"].isna().all()

    summary = condition_summary(df)
    assert summary.columns.tolist() == ["Condition", "Samples", "Mean Quantified proteins"]
    assert summary.set_index("Condition")["Mean Quantified proteins"].to_dict() == {"A": 1.5, "B": 3.0}
