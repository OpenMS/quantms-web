"""The source of the Viewer's charge state and missed cleavage component: the
IDFilter idXMLs (target PSMs at the run's PSM FDR) turned into the columns the
library's ``src.viz.identification.charge_cleavage`` reads (``file``,
``sequence``, ``charge``, ``missed_cleavages``).

Missed cleavages are counted for the enzyme each file was searched with (its
search parameters), so a run digested with another protease counts its own
cleavage sites, and left empty for a search without an enzyme.
Streamlit-free, so tests call it.
"""
from pathlib import Path

import pandas as pd
import pyopenms as poms

from src.common.qc_plots import _is_decoy
from src.common.results_helpers import extract_filename_from_idxml, load_idxml

DEFAULT_ENZYME = "Trypsin/P"
# Enzymes without cleavage sites: a missed cleavage means nothing there
NO_CLEAVAGE_SITES = {"unspecific cleavage", "no cleavage", "unspecific", "none", ""}

PSM_COLUMNS = ["file", "sequence", "charge", "missed_cleavages"]


def searched_enzyme(protein_ids) -> str:
    """The enzyme the file was searched with, from its search parameters."""
    for prot in protein_ids:
        try:
            name = prot.getSearchParameters().digestion_enzyme.getName()
        except (AttributeError, RuntimeError):
            continue
        # A file written without search parameters names "unknown_enzyme"
        if name and name != "unknown_enzyme":
            return name
    return DEFAULT_ENZYME


def _digestion(enzyme: str):
    """A ProteaseDigestion for ``enzyme``, None for an enzyme without cleavage sites."""
    if enzyme.lower() in NO_CLEAVAGE_SITES:
        return None
    digestion = poms.ProteaseDigestion()
    try:
        digestion.setEnzyme(enzyme)
    except RuntimeError:
        digestion.setEnzyme(DEFAULT_ENZYME)
    return digestion


def read_psms(idxml_path: Path, mzml_names=None) -> tuple[pd.DataFrame, str]:
    """The target PSMs (top hit of each spectrum) of one idXML: file, sequence,
    charge and missed cleavages, and the enzyme the file was searched with.
    Missed cleavages are empty for a search without an enzyme."""
    proteins, peptides = load_idxml(idxml_path)
    enzyme = searched_enzyme(proteins)
    digestion = _digestion(enzyme)
    file = extract_filename_from_idxml(Path(idxml_path), mzml_names)
    rows = []
    for pep in peptides:
        hits = pep.getHits()
        if not hits or _is_decoy(hits[0]):
            continue
        hit = hits[0]
        unmodified = hit.getSequence().toUnmodifiedString()
        rows.append({
            "file": file,
            "sequence": unmodified,
            "charge": hit.getCharge(),
            "missed_cleavages": digestion.countInternalCleavageSites(unmodified) if digestion else None,
        })
    return pd.DataFrame(rows, columns=PSM_COLUMNS), enzyme
