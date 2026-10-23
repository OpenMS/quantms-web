"""Common contaminants (cRAP) for the search database.

The workflow appends the cRAP proteins to the user's FASTA with ``CON__`` before
their accessions, so every protein and PSM table can flag them, and the
protein-level results and downstream analysis leave them out by default.
"""
from pathlib import Path

import pandas as pd
from pyopenms import FASTAFile

# Accession prefix of the contaminants appended to the database (as MaxQuant and ProSE)
CONTAMINANT_PREFIX = "CON__"

# The common Repository of Adventitious Proteins (cRAP, thegpm.org/crap), 116 entries
CRAP_FASTA = Path("example-data", "db", "contaminants", "crap.fasta")

# Separators between the accessions of a protein group or of a PSM's proteins
_ACCESSION_SEPARATORS = (";", ",", "/")


def write_contaminants(out_path, crap_path=CRAP_FASTA) -> int:
    """Write the cRAP entries with ``CON__`` accessions to ``out_path``.

    Entries that already carry the prefix keep it. Returns the number of entries.
    """
    entries = []
    FASTAFile().load(str(crap_path), entries)
    for entry in entries:
        if not entry.identifier.startswith(CONTAMINANT_PREFIX):
            entry.identifier = CONTAMINANT_PREFIX + entry.identifier
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    FASTAFile().store(str(out_path), entries)
    return len(entries)


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
