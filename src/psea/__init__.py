"""PTM enrichment: which protein modifications are enriched among the proteins
that change between two groups.

The proteins of the group comparison, ranked by their log2 fold change, go to
protein set enrichment analysis (PSEA) of the PEIMAN2 R package (``psea.R``),
which tests each modification UniProt annotates for the organism. The workflow
runs it at the end of a run, on the comparison with the recommended settings
(``WorkflowTest.ptm_enrichment``); the Viewer shows the result and runs it
again on the comparison as edited.

Streamlit-free but ``run_scores``, which reads the comparison the way the
Viewer does.
"""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).with_name("psea.R")
RESULTS = "ptm_enrichment"
SCORES_FILE = "protein_scores.csv"
RESULTS_FILE = "psea_results.csv"

# Organisms exactly as PEIMAN2's bundled table names them (its OS column; psea.R
# also resolves other spellings), with
# the species code of their UniProt entry names (``ALBU_HUMAN``)
_KNOWN = [
    ("HUMAN", "Homo sapiens (Human)"),
    ("MOUSE", "Mus musculus (Mouse)"),
    ("RAT", "Rattus norvegicus (Rat)"),
    ("BOVIN", "Bos taurus (Bovine)"),
    ("PIG", "Sus scrofa (Pig)"),
    ("CANLF", "Canis lupus familiaris (Dog) (Canis familiaris)"),
    ("CHICK", "Gallus gallus (Chicken)"),
    ("XENLA", "Xenopus laevis (African clawed frog)"),
    ("DANRE", "Danio rerio (Zebrafish) (Brachydanio rerio)"),
    ("DROME", "Drosophila melanogaster (Fruit fly)"),
    ("CAEEL", "Caenorhabditis elegans"),
    ("ARATH", "Arabidopsis thaliana (Mouse-ear cress)"),
    ("ORYSJ", "Oryza sativa subsp japonica (Rice)."),
    ("YEAST", "Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)"),
    ("SCHPO", "Schizosaccharomyces pombe (strain 972 / ATCC 24843) (Fission yeast)"),
    ("PLAF7", "Plasmodium falciparum (isolate 3D7)"),
    ("ECOLI", "Escherichia coli (strain K12)"),
    ("BACSU", "Bacillus subtilis (strain 168)"),
    ("MYCTU", "Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv)"),
]
ORGANISMS = [name for _, name in _KNOWN]
DEFAULT_ORGANISM = ORGANISMS[0]
DEFAULT_PERMUTATIONS = 1000
# The choices of the Organism menu besides the organisms
AUTO = "Detect automatically"
OTHER = "Other (type its UniProt name)"
CHOICES = [AUTO, *ORGANISMS, OTHER]
_MNEMONICS = dict(_KNOWN)
_HEADERS_READ = 5000  # FASTA entries read to detect the organism


def _scientific(name: str) -> str:
    """``homo sapiens`` of ``Homo sapiens (Human)``: lower case, without
    full stops (UniProt writes ``subsp.``, PEIMAN2's table ``subsp``)."""
    return re.split(r"\s*\(", str(name), maxsplit=1)[0].replace(".", "").strip().lower()


def organism_named(name) -> str | None:
    """The UniProt organism name of ``name``: a listed organism when its
    scientific name matches (``Homo sapiens`` → ``Homo sapiens (Human)``),
    else ``name`` as it is (None: empty)."""
    name = str(name or "").strip()
    if not name:
        return None
    wanted = _scientific(name)
    for organism in ORGANISMS:
        if _scientific(organism) == wanted:
            return organism
    return name


def settings(params: dict, database=None) -> tuple[bool, str | None, int]:
    """(switched on, organism, permutations) of a run's or the next run's
    settings. The organism is the one chosen in the Organism menu; with
    *Detect automatically* (and in settings saved before the menu had it,
    without a choice) the one ``database`` (the protein FASTA) names, None
    when it names none."""
    choice = params.get("ptm-organism") or AUTO
    if choice == OTHER:
        organism = organism_named(params.get("ptm-organism-other"))
    elif choice == AUTO:
        organism = detect(database)
    else:
        organism = choice
    try:
        permutations = max(1, int(params.get("ptm-permutations", DEFAULT_PERMUTATIONS)))
    except (TypeError, ValueError):
        permutations = DEFAULT_PERMUTATIONS
    return bool(params.get("ptm-enrichment", False)), organism, permutations


def guess_organism(names) -> str | None:
    """The organism UniProt entry names (``ALBU_HUMAN``, ``sp|P02768|ALBU_HUMAN``)
    or a database file name point to, by majority (None: no hint)."""
    counts: dict = {}
    for name in names:
        upper = str(name).upper()
        for mnemonic, organism in _MNEMONICS.items():
            if re.search(rf"(?:^|[_|]){mnemonic}(?:$|[^A-Z0-9])", upper) or upper.startswith(mnemonic + "_"):
                counts[organism] = counts.get(organism, 0) + 1
    return max(counts, key=counts.get) if counts else None


def fasta_organism(path) -> str | None:
    """The organism most entries of a FASTA file belong to: the ``OS=`` field
    of UniProt headers, else the species code of their entry names (None:
    neither, or no readable file). Reads the first entries only."""
    try:
        handle = open(path, encoding="utf-8", errors="replace")
    except (OSError, TypeError):
        return None
    species: dict = {}
    entries = []
    with handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if header.upper().startswith(("DECOY", "REV_")):
                continue
            entries.append(header.split()[0] if header else "")
            found = re.search(r"\bOS=(.+?)(?:\s+[A-Z]{2}=|$)", header)
            if found:
                species[found.group(1).strip()] = species.get(found.group(1).strip(), 0) + 1
            if len(entries) >= _HEADERS_READ:
                break
    if species:
        return organism_named(max(species, key=species.get))
    return guess_organism(entries)


def detect(database) -> str | None:
    """The organism of a protein database: from its entries, else its file name."""
    if not database:
        return None
    path = Path(str(database))
    try:
        stamp = path.stat().st_mtime_ns
    except OSError:
        stamp = None
    return _detect(str(path), stamp)


@lru_cache(maxsize=32)
def _detect(path: str, stamp) -> str | None:
    found = fasta_organism(path) if stamp is not None else None
    return found or guess_organism([Path(path).name])


def uniprot_accession(identifier) -> str:
    """The UniProt accession (``P02768``) in ``sp|P02768|ALBU_HUMAN``; other identifiers as they are."""
    parts = str(identifier).strip().split("|")
    return parts[1].strip() if len(parts) >= 3 and parts[0] in ("sp", "tr") else str(identifier).strip()


def scores(statistics_df: pd.DataFrame, id_col: str) -> pd.DataFrame:
    """The ranked list PSEA reads: ``UniProtAC, Score`` (the UniProt accession
    of the lead protein of each compared group and its log2 fold change;
    contaminants and decoys left out, a protein listed twice kept with its
    largest change)."""
    from src.common.contaminants import is_contaminant
    from src.common.differential import lead_accession

    if statistics_df is None or statistics_df.empty or "log2FC" not in statistics_df:
        return pd.DataFrame(columns=["UniProtAC", "Score"])
    groups = statistics_df[id_col].astype(str)
    keep = ~groups.map(is_contaminant)
    leads = groups.map(lead_accession)
    keep &= ~leads.str.upper().str.startswith(("DECOY", "REV_"))
    table = pd.DataFrame({
        "UniProtAC": leads.map(uniprot_accession),
        "Score": pd.to_numeric(statistics_df["log2FC"], errors="coerce"),
    })[keep]
    table = table[table["UniProtAC"].astype(bool)]
    table = table[table["Score"].map(lambda v: pd.notna(v) and abs(v) != float("inf"))]
    if table.empty:
        return pd.DataFrame(columns=["UniProtAC", "Score"])
    table = table.loc[table["Score"].abs().groupby(table["UniProtAC"]).idxmax()]
    return table.sort_values("Score", ascending=False).reset_index(drop=True)


def command(out_dir: Path, scores_csv: Path, organism: str, permutations: int) -> list[str]:
    return [rscript() or "Rscript", str(SCRIPT), str(out_dir), str(scores_csv), organism, str(int(permutations))]


def rscript() -> str | None:
    return shutil.which("Rscript")


@lru_cache(maxsize=1)
def missing() -> str | None:
    """Why PTM enrichment cannot run here (None: it can)."""
    if rscript() is None:
        return "R is not installed here."
    try:
        found = subprocess.run(
            [rscript(), "-e", "quit(status = if (requireNamespace('PEIMAN2', quietly = TRUE)) 0 else 1)"],
            capture_output=True, timeout=60, check=False,
        ).returncode == 0
    except (OSError, subprocess.SubprocessError):
        found = False
    return None if found else "The R package PEIMAN2 is not installed here."


def digest(table: pd.DataFrame, organism: str, permutations: int) -> str:
    return hashlib.sha1(
        repr((table.to_csv(index=False), organism, int(permutations))).encode()
    ).hexdigest()[:12]


def write_scores(table: pd.DataFrame, folder: Path) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / SCORES_FILE
    table.to_csv(path, index=False)
    return path


def compute(table: pd.DataFrame, organism: str, permutations: int, folder: Path) -> tuple[Path | None, str | None]:
    """Run PSEA on ``table`` into ``folder`` (kept: a second call with the same
    input reads it). Returns (the result folder, None) or (None, why not)."""
    folder = Path(folder) / digest(table, organism, permutations)
    if (folder / RESULTS_FILE).is_file():
        return folder, None
    why = missing()
    if why:
        return None, why
    if table.empty:
        return None, "No compared protein has a fold change."
    scores_csv = write_scores(table, folder)
    done = subprocess.run(command(folder, scores_csv, organism, permutations),
                          capture_output=True, text=True, check=False)
    if done.returncode != 0 or not (folder / RESULTS_FILE).is_file():
        lines = [l for l in (done.stderr or done.stdout or "").splitlines() if l.strip()]
        return None, lines[-1] if lines else f"Rscript stopped with exit code {done.returncode}."
    return folder, None


def read_results(folder: Path | None) -> pd.DataFrame | None:
    path = Path(folder, RESULTS_FILE) if folder else None
    if path is None or not path.is_file():
        return None
    try:
        return pd.read_csv(path)
    except (pd.errors.ParserError, pd.errors.EmptyDataError, OSError):
        return None


def run_results(results_dir: Path) -> Path | None:
    """The folder of the result the workflow computed for a run (None: none)."""
    folder = Path(results_dir, RESULTS)
    return folder if (folder / RESULTS_FILE).is_file() else None


def run_scores(workspace: Path, run_id: str, params: dict) -> pd.DataFrame:
    """The ranked list of a finished run: its group comparison with the
    recommended settings, read as the Viewer reads it (the workflow computed
    it in ``precompute`` just before)."""
    import streamlit as st

    from src.common.guide import ensure_group_comparison
    from src.common.results_helpers import get_abundance_data, get_id_column
    from src.viewer.scope import experiment

    previous = st.session_state.get("workspace")
    st.session_state["workspace"] = Path(workspace)
    try:
        with experiment(run_id, 0, single=True):
            statistics_df = ensure_group_comparison(params, show_errors=False)
            abundance = get_abundance_data(Path(workspace))
            if statistics_df is None or abundance is None:
                return scores(None, "")
            return scores(statistics_df, get_id_column(Path(workspace), abundance[0]))
    finally:
        if previous is None:
            st.session_state.pop("workspace", None)
        else:
            st.session_state["workspace"] = previous
