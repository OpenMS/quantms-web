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
    assert not on and organism == psea.DEFAULT_ORGANISM and permutations == psea.DEFAULT_PERMUTATIONS
    on, organism, permutations = psea.settings({"ptm-enrichment": True, "ptm-organism": psea.ORGANISMS[1],
                                                "ptm-permutations": 200})
    assert on and organism == psea.ORGANISMS[1] and permutations == 200


def test_command_runs_the_script():
    cmd = psea.command("/out", "/out/scores.csv", "Homo sapiens (Human)", 500)
    assert cmd[1] == str(psea.SCRIPT) and cmd[2:] == ["/out", "/out/scores.csv", "Homo sapiens (Human)", "500"]
