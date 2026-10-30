"""Common contaminants (cRAP) for the search database.

The workflow appends the cRAP proteins to the user's FASTA with ``CON__`` before
their accessions, so every protein and PSM table can flag them, and the
protein-level results and downstream analysis leave them out by default.
"""
import re
from pathlib import Path

import pandas as pd
from pyopenms import FASTAFile

# Accession prefix of the contaminants appended to the database (as MaxQuant and ProSE)
CONTAMINANT_PREFIX = "CON__"

# The common Repository of Adventitious Proteins (cRAP, thegpm.org/crap), 116 entries
CRAP_FASTA = Path("example-data", "db", "contaminants", "crap.fasta")

# Separators between the accessions of a protein group or of a PSM's proteins
_ACCESSION_SEPARATORS = (";", ",", "/")


def write_contaminants(out_path, crap_path=CRAP_FASTA, targets=None) -> tuple[int, int]:
    """Write the cRAP entries with ``CON__`` accessions to ``out_path``.

    Entries that already carry the prefix keep it. A cRAP protein that is in the
    ``targets`` FASTA, by sequence or by UniProt entry name (cRAP's sequences
    can differ from the current UniProt ones), is left out: the user searches for it on purpose
    (e.g. the UPS spike-in standards, which cRAP lists too), and a CON__ twin
    would mark its protein group as a contaminant and drop it from the results.
    The file is rewritten only when its content changes, so a decoy database
    built from it stays valid. Returns ``(written, already_in_targets)``.
    """
    entries = []
    FASTAFile().load(str(crap_path), entries)
    target_sequences, target_names = set(), set()
    if targets is not None:
        target_entries = []
        FASTAFile().load(str(targets), target_entries)
        target_sequences = {e.sequence for e in target_entries}
        target_names = {_entry_name(e.identifier) for e in target_entries} - {""}
    kept = [
        e for e in entries
        if e.sequence not in target_sequences and _entry_name(e.identifier) not in target_names
    ]
    for entry in kept:
        if not entry.identifier.startswith(CONTAMINANT_PREFIX):
            entry.identifier = CONTAMINANT_PREFIX + entry.identifier
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".tmp")
    FASTAFile().store(str(tmp), kept)
    if out_path.exists() and out_path.read_bytes() == tmp.read_bytes():
        tmp.unlink()
    else:
        tmp.replace(out_path)
    return len(kept), len(entries) - len(kept)


def _entry_name(identifier: str) -> str:
    """The UniProt entry name (``ALBU_HUMAN``) in ``sp|P02768|ALBU_HUMAN`` or cRAP's ``sp|ALBU_HUMAN|``."""
    for part in identifier.removeprefix(CONTAMINANT_PREFIX).split("|"):
        if re.fullmatch(r"[A-Z0-9]+_[A-Z0-9]+", part):
            return part
    return ""


def append_fasta(parts, out_path) -> None:
    """Write the entries of the FASTA files ``parts``, in order, to ``out_path``."""
    entries = []
    for part in parts:
        part_entries = []
        FASTAFile().load(str(part), part_entries)
        entries.extend(part_entries)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    FASTAFile().store(str(out_path), entries)


def is_contaminant(accessions) -> bool:
    """Whether a protein, protein group or PSM maps to a contaminant.

    ``accessions`` is one accession or several joined by ``;``, ``,`` or ``/``
    (ProteinQuantifier groups). A group counts as a contaminant when any of its
    members is one, as in MaxQuant.
    """
    if not isinstance(accessions, str) or not accessions:
        return False
    for sep in _ACCESSION_SEPARATORS[1:]:
        accessions = accessions.replace(sep, _ACCESSION_SEPARATORS[0])
    return any(a.strip().startswith(CONTAMINANT_PREFIX) for a in accessions.split(_ACCESSION_SEPARATORS[0]))


def contaminant_mask(accessions: pd.Series) -> pd.Series:
    """Boolean Series: True where the row's accessions include a contaminant."""
    return accessions.map(is_contaminant).astype(bool)


def drop_contaminants(df: pd.DataFrame, column: str) -> pd.DataFrame:
    """``df`` without the rows whose ``column`` names a contaminant."""
    if df is None or column not in df.columns:
        return df
    return df[~contaminant_mask(df[column])]
