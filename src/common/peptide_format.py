"""Peptide sequences with modifications, written for reading or for export.

OpenMS writes a modified peptide as ``.(TMT6plex)PEPM(Oxidation)K``. Two other
forms are offered here:

- labels: each modified residue carries a short label, ``tmt-PEPM(ox)K``, for
  scanning a table by eye;
- ProForma 2.0 (HUPO-PSI), ``[TMT6plex]-PEPM[Oxidation]K``, the standard form
  other tools read, and what exports use.

Both are built with pyOpenMS from the OpenMS string, so any modification
OpenMS knows is handled; one without a Unimod name falls back to its mass shift.
It makes no Streamlit calls.
"""
from functools import lru_cache

from pyopenms import AASequence

# Short labels for common Unimod modifications; any other one shows its full name
SHORT_LABELS = {
    "Oxidation": "ox",
    "Phospho": "ph",
    "Carbamidomethyl": "cam",
    "Acetyl": "ac",
    "Deamidated": "de",
    "Methyl": "me",
    "Dimethyl": "me2",
    "Trimethyl": "me3",
    "GG": "gg",
    "Gln->pyro-Glu": "pyro",
    "Glu->pyro-Glu": "pyro",
    "Ammonia-loss": "-nh3",
    "Carbamyl": "carb",
    "Nitro": "no2",
    "Amidated": "am",
    "TMT6plex": "tmt",
    "TMTpro": "tmtpro",
    "iTRAQ4plex": "itraq",
    "iTRAQ8plex": "itraq8",
}

DISPLAY_MODES = ("Labels", "ProForma", "OpenMS")


def _name(mod) -> str:
    """Unimod name of a modification, or its signed mass shift when it has none."""
    name = mod.getId()
    if name:
        return name
    return f"{mod.getDiffMonoMass():+.4f}"


def _modified_parts(sequence: AASequence):
    """(N-term mod or None, [(residue letter, mod or None)], C-term mod or None)."""
    n_term = sequence.getNTerminalModification() if sequence.hasNTerminalModification() else None
    c_term = sequence.getCTerminalModification() if sequence.hasCTerminalModification() else None
    residues = []
    for i in range(sequence.size()):
        residue = sequence.getResidue(i)
        mod = residue.getModification() if residue.isModified() else None
        residues.append((residue.getOneLetterCode(), mod))
    return n_term, residues, c_term


@lru_cache(maxsize=65536)
def to_proforma(openms_sequence: str) -> str:
    """ProForma 2.0 form of an OpenMS peptide string; unparsable strings come back as they are."""
    try:
        sequence = AASequence.fromString(openms_sequence)
    except Exception:  # not a peptide pyOpenMS can read: show it unchanged
        return openms_sequence
    n_term, residues, c_term = _modified_parts(sequence)
    text = "".join(f"{aa}[{_name(mod)}]" if mod else aa for aa, mod in residues)
    if n_term:
        text = f"[{_name(n_term)}]-{text}"
    if c_term:
        text = f"{text}-[{_name(c_term)}]"
    return text


@lru_cache(maxsize=65536)
def to_labels(openms_sequence: str) -> str:
    """Residues followed by short modification labels, e.g. ``tmt-PEPM(ox)K``."""
    try:
        sequence = AASequence.fromString(openms_sequence)
    except Exception:  # not a peptide pyOpenMS can read: show it unchanged
        return openms_sequence

    def label(mod) -> str:
        name = _name(mod)
        return SHORT_LABELS.get(name, name)

    n_term, residues, c_term = _modified_parts(sequence)
    text = "".join(f"{aa}({label(mod)})" if mod else aa for aa, mod in residues)
    if n_term:
        text = f"{label(n_term)}-{text}"
    if c_term:
        text = f"{text}-{label(c_term)}"
    return text


def format_sequence(openms_sequence: str, mode: str) -> str:
    """``openms_sequence`` in one of DISPLAY_MODES."""
    if mode == "ProForma":
        return to_proforma(openms_sequence)
    if mode == "Labels":
        return to_labels(openms_sequence)
    return openms_sequence
