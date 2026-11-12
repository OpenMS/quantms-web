"""Peak memory of the ProSE database search of a run (Streamlit-free).

ProSE builds a fragment index of the database (targets and decoys) and keeps
it while it searches its files one after the other; the file being searched is
read whole. Several ProSE processes run side by side (``search_batches``), each
with an index of its own. The estimate per process:

    BASE
    + index:   INDEX_PER_1000_PROTEINS x proteins x SPECIFICITY x MOD^(variable mods - 1)
               non-specific: PER_MILLION_RESIDUES x residues x MOD^...
               chunked (database:chunk_size below the proteins): the index of one chunk,
               which costs CHUNKED times the unchunked index per protein, at most the whole
    + 2 x the largest mzML of the process + the FASTA
    + growth: each further file of the process adds GROWTH (levels off after the 5th)
    + THREAD x every thread above 4

The constants are measured peaks of ProSE (a 64-bit WebAssembly build, whose
allocator was set up to match the native one) over tryptic, semi-tryptic and
non-specific searches of a bacterial and the human proteome (1 to 12 files of
75 MB to 730 MB, 4 to 16 threads), at or above every measured peak. They are a
calibration, not a bound: the estimate is meant to warn, never to refuse. A protein-terminal variable
modification does not enlarge the index (only one peptide per protein carries
it), so it is not counted.
"""

from functools import lru_cache
from pathlib import Path

GB = 1024 ** 3
MB = 1024 ** 2

BASE = 0.4 * GB
INDEX_PER_1000_PROTEINS = 0.025 * GB
SPECIFICITY = {"full": 1.0, "semi": 15.0}
CHUNKED = {"full": 32 / 15, "semi": 32.0}
PER_MILLION_RESIDUES = 9.6 * GB
PER_MILLION_CHUNK_RESIDUES = 24.5 * GB
MOD = 2.0
GROWTH_SMALL, GROWTH_LARGE = 0.06 * GB, 0.15 * GB
SMALL_FILE, LARGE_FILE = 75 * MB, 733 * MB
THREAD = 0.035 * GB
RESIDUES_PER_PROTEIN = 555  # where a database's are unknown (human's)


@lru_cache(maxsize=16)
def _fasta_counts(path: str, size: int, mtime: float) -> tuple[int, int]:
    proteins = residues = 0
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                proteins += 1
            else:
                residues += len(line.strip())
    return proteins, residues


def fasta_counts(path) -> tuple[int, int]:
    """(proteins, residues) of a FASTA file; read once per size and time."""
    p = Path(path)
    st = p.stat()
    return _fasta_counts(str(p), st.st_size, st.st_mtime)


def index_modifications(variable) -> int:
    """Variable modifications that enlarge the index (not protein-terminal ones)."""
    if isinstance(variable, str):
        variable = variable.split("\n")
    names = [str(v).strip() for v in variable or [] if str(v).strip()]
    return sum(1 for n in names if not (n.endswith("(Protein N-term)") or n.endswith("(Protein C-term)")))


def specificity(values: dict) -> str:
    if values.get("Search:enzyme") == "unspecific cleavage":
        return "none"
    spec = str(values.get("Search:peptide:enzyme_specificity", "full"))
    return spec if spec in ("full", "semi", "none") else "full"


def index_bytes(proteins: int, residues: int, spec: str, mods: int, chunk: int = 0) -> float:
    """The fragment index of one ProSE process."""
    factor = MOD ** max(0, mods - 1)
    chunked = 0 < chunk < proteins
    if spec == "none":
        residues = residues or proteins * RESIDUES_PER_PROTEIN
        whole = PER_MILLION_RESIDUES * factor * residues / 1e6
        if not chunked:
            return whole
        chunk_residues = chunk * residues / max(1, proteins)
        return min(PER_MILLION_CHUNK_RESIDUES * factor * chunk_residues / 1e6, whole)
    per_protein = INDEX_PER_1000_PROTEINS / 1000 * factor
    whole = per_protein * proteins * SPECIFICITY[spec]
    return min(per_protein * chunk * CHUNKED[spec], whole) if chunked else whole


def growth(files: int, mean_bytes: float) -> float:
    """Heap growth of the files after the first of one process."""
    at = min(1.0, max(0.0, (mean_bytes - SMALL_FILE) / (LARGE_FILE - SMALL_FILE)))
    g = GROWTH_SMALL + (GROWTH_LARGE - GROWTH_SMALL) * at
    return g * min(max(0, files - 1), 4) + g / 4 * max(0, files - 5)


def search_memory(file_sizes: list, batches: list, fasta_bytes: int, proteins: int,
                  residues: int, values: dict, threads: int) -> dict:
    """Peak memory of the search: every process of ``batches`` at once.

    Returns ``{"bytes", "index", "inputs", "cause"}``; ``cause`` is what
    dominates, ``"database"`` (the index) or ``"files"``.
    """
    spec = specificity(values)
    mods = index_modifications(values.get("Search:modifications:variable", ""))
    try:
        chunk = int(values.get("Search:database:chunk_size", 0) or 0)
    except (TypeError, ValueError):
        chunk = 0
    index = index_bytes(proteins, residues, spec, mods, chunk)
    per_process_threads = max(1, int(threads) // max(1, len(batches)))
    total = inputs = 0.0
    for batch in batches:
        sizes = [file_sizes[i] for i in batch] or [0]
        files_part = 2 * max(sizes) + fasta_bytes + growth(len(sizes), sum(sizes) / len(sizes))
        inputs += files_part
        total += BASE + index + files_part + THREAD * max(0, per_process_threads - 4)
    return {
        "bytes": int(total),
        "index": int(index * len(batches)),
        "inputs": int(inputs),
        "cause": "database" if index * len(batches) >= inputs else "files",
    }
