"""PTM enrichment: the ranked protein list handed to the enrichment, and its settings."""
import pandas as pd

from src import psea


def test_scores_rank_lead_accessions_by_fold_change():
    stats = pd.DataFrame({
        "ProteinName": ["sp|P1|A_HUMAN;sp|P9|Z_HUMAN", "sp|P2|B_HUMAN", "CON__sp|P3|C_HUMAN", "DECOY_sp|P4|D",
                        "sp|P5|E_HUMAN", "sp|P1|A_HUMAN"],
        "log2FC": [1.5, -2.0, 3.0, 4.0, float("nan"), -0.5],
    })
    table = psea.scores(stats, "ProteinName")
    assert list(table.columns) == ["UniProtAC", "Score"]
    assert table["UniProtAC"].tolist() == ["P1", "P2"]
    assert table["Score"].tolist() == [1.5, -2.0]


def test_settings_default_off():
    on, organism, permutations = psea.settings({})
    assert not on and organism is None and permutations == psea.DEFAULT_PERMUTATIONS
    on, organism, permutations = psea.settings({"ptm-enrichment": True, "ptm-organism": psea.ORGANISMS[1],
                                                "ptm-permutations": 200})
    assert on and organism == psea.ORGANISMS[1] and permutations == 200


def test_command_runs_the_script():
    cmd = psea.command("/out", "/out/scores.csv", "Homo sapiens (Human)", 500)
    assert cmd[1] == str(psea.SCRIPT) and cmd[2:] == ["/out", "/out/scores.csv", "Homo sapiens (Human)", "500"]


def _fasta(tmp_path, name, headers):
    path = tmp_path / name
    path.write_text("".join(f">{h}\nPEPTIDEK\n" for h in headers))
    return path


def test_detect_reads_uniprot_headers_not_the_file_name(tmp_path):
    db = _fasta(tmp_path, "database.fasta", [
        "sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2",
        "sp|P07724|ALBU_MOUSE Albumin OS=Mus musculus OX=10090 GN=Alb PE=1 SV=3",
        "sp|P00001|X_MOUSE Something OS=Mus musculus OX=10090",
        "DECOY_sp|P00002|Y_HUMAN decoy OS=Homo sapiens OX=9606",
    ])
    assert psea.detect(db) == "Mus musculus (Mouse)"


def test_detect_yeast_strain_and_unlisted_organism(tmp_path):
    yeast = _fasta(tmp_path, "y.fasta", [
        "sp|P00330|ADH1_YEAST Alcohol dehydrogenase 1 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=ADH1"])
    assert psea.detect(yeast) == psea.ORGANISMS[psea.ORGANISMS.index(
        "Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)")]
    other = _fasta(tmp_path, "o.fasta", ["tr|A0A000|A0A000_9ZZZZ Protein OS=Ovis aries OX=9940"])
    assert psea.detect(other) == "Ovis aries"


def test_detect_falls_back_to_entry_names_then_file_name(tmp_path):
    names = _fasta(tmp_path, "db.fasta", ["sp|P02768|ALBU_HUMAN", "sp|P00001|A_YEAST", "sp|P00002|B_YEAST"])
    assert psea.detect(names) == psea.ORGANISMS[psea.ORGANISMS.index(
        "Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)")]
    plain = _fasta(tmp_path, "uniprot_rat_2024.fasta", ["P1 some protein", "P2 another"])
    assert psea.detect(plain) == "Rattus norvegicus (Rat)"
    assert psea.detect(_fasta(tmp_path, "db2.fasta", ["P1 some protein"])) is None
    assert psea.detect(None) is None


def test_settings_follow_the_organism_menu(tmp_path):
    db = _fasta(tmp_path, "db.fasta", ["sp|P1|A_PIG x OS=Sus scrofa OX=9823"])
    assert psea.settings({"ptm-organism": psea.AUTO}, db)[1] == "Sus scrofa (Pig)"
    assert psea.settings({}, db)[1] == "Sus scrofa (Pig)"
    assert psea.settings({"ptm-organism": psea.ORGANISMS[2]}, db)[1] == psea.ORGANISMS[2]
    assert psea.settings({"ptm-organism": psea.OTHER, "ptm-organism-other": "mus musculus"})[1] == \
        "Mus musculus (Mouse)"
    assert psea.settings({"ptm-organism": psea.OTHER, "ptm-organism-other": "Ovis aries (Sheep)"})[1] == \
        "Ovis aries (Sheep)"
    assert psea.settings({"ptm-organism": psea.OTHER})[1] is None


def test_guess_organism_from_protein_names():
    assert psea.guess_organism(["sp|P1|A_HUMAN;sp|P2|B_HUMAN", "C_MOUSE"]) == "Homo sapiens (Human)"
    assert psea.guess_organism(["RATIO_1", "PIGMENT"]) is None
