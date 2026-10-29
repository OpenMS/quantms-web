"""Protein-level identification tables: proteins, their peptides and PSMs, and sequence coverage.

The tables are built from the FDR-filtered PSMs of every run and the target
FASTA database, so they exist for LFQ and TMT alike. Decoy hits are left out.
"""
import html
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pyopenms import FASTAFile

from src.common.results_helpers import extract_filename_from_idxml, extract_scan_from_ref, load_idxml

DECOY_PREFIXES = ("DECOY_", "decoy_", "rev_", "REV_", "XXX_", "reversed_")
Q_VALUE_KEYS = ("MS:1001491", "q-value")

PSM_COLUMNS = ["file", "scan", "rt", "mz", "charge", "peptide", "sequence", "score", "q_value", "proteins"]


def _is_decoy_accession(accession: str) -> bool:
    return accession.startswith(DECOY_PREFIXES)


def read_fasta(path) -> dict[str, tuple[str, str]]:
    """Accession -> (description, sequence) of the target entries of a FASTA file."""
    entries = []
    FASTAFile().load(str(path), entries)
    return {
        e.identifier: (e.description, e.sequence)
        for e in entries
        if not _is_decoy_accession(e.identifier)
    }


def gene_name(description: str) -> str:
    """The gene from a UniProt description (``GN=``), or an empty string."""
    match = re.search(r"\bGN=(\S+)", description or "")
    return match.group(1) if match else ""


def protein_name(description: str) -> str:
    """A UniProt description without its ``OS=``/``OX=``/``GN=``/``PE=``/``SV=`` fields."""
    return re.split(r"\s+[A-Z]{2}=", description or "", maxsplit=1)[0]


def collect_psms(idxml_paths, mzml_names=None) -> pd.DataFrame:
    """Top hits of all runs, one row per target PSM, with their target proteins."""
    records = []
    for path in idxml_paths:
        file = extract_filename_from_idxml(Path(path), mzml_names)
        for pep in load_idxml(path)[1]:
            hits = pep.getHits()
            if not hits:
                continue
            hit = hits[0]
            proteins = sorted({
                ev.getProteinAccession() for ev in hit.getPeptideEvidences()
                if not _is_decoy_accession(ev.getProteinAccession())
            })
            if not proteins:
                continue  # a decoy hit
            ref = pep.getMetaValue("spectrum_reference") if pep.metaValueExists("spectrum_reference") else ""
            if isinstance(ref, bytes):
                ref = ref.decode()
            q_value = np.nan
            for key in Q_VALUE_KEYS:
                if hit.metaValueExists(key):
                    q_value = float(hit.getMetaValue(key))
                    break
            seq = hit.getSequence()
            records.append({
                "file": file,
                "scan": extract_scan_from_ref(ref),
                "rt": pep.getRT(),
                "mz": pep.getMZ(),
                "charge": hit.getCharge(),
                "peptide": seq.toString(),
                "sequence": seq.toUnmodifiedString(),
                "score": hit.getScore(),
                "q_value": q_value,
                "proteins": proteins,
            })
    return pd.DataFrame(records, columns=PSM_COLUMNS)


def coverage_mask(protein_sequence: str, peptides) -> np.ndarray:
    """Per residue, whether some peptide covers it (every occurrence of each peptide counts)."""
    mask = np.zeros(len(protein_sequence), dtype=bool)
    for peptide in set(peptides):
        if not peptide:
            continue
        start = protein_sequence.find(peptide)
        while start != -1:
            mask[start:start + len(peptide)] = True
            start = protein_sequence.find(peptide, start + 1)
    return mask


def build_tables(psms: pd.DataFrame, fasta: dict[str, tuple[str, str]], files: list[str]):
    """Proteins, peptide-to-protein and PSM-to-protein tables.

    Peptides and PSMs are in long format (one row per protein they map to), so
    selecting a protein filters them by its accession.

    Returns:
        (proteins, peptides, psms_long, coverage) where ``coverage`` maps an
        accession to its residue mask.
    """
    long = psms.explode("proteins").rename(columns={"proteins": "accession"})
    long["shared"] = psms["proteins"].str.len().reindex(long.index).gt(1)

    rows, coverage = [], {}
    for accession, group in long.groupby("accession", sort=False):
        description, sequence = fasta.get(accession, ("", ""))
        mask = coverage_mask(sequence, group["sequence"])
        coverage[accession] = mask
        per_file = group.groupby("file").size()
        row = {
            "accession": accession,
            "gene": gene_name(description),
            "description": protein_name(description),
            "coverage (%)": round(100 * mask.mean(), 1) if len(mask) else np.nan,
            "peptides": group["sequence"].nunique(),
            "unique peptides": group.loc[~group["shared"], "sequence"].nunique(),
            "PSMs": len(group),
            "found in": f"{per_file.size} of {len(files)}",
            "length": len(sequence) or np.nan,
        }
        for file in files:
            row[f"PSMs {Path(file).stem}"] = int(per_file.get(file, 0))
        rows.append(row)
    proteins = pd.DataFrame(rows)
    if not proteins.empty:
        proteins = proteins.sort_values(["PSMs", "peptides"], ascending=False, ignore_index=True)
        # the description last: it is the widest column, and the one to truncate
        proteins = proteins[[c for c in proteins.columns if c != "description"] + ["description"]]
        if not proteins["gene"].astype(bool).any():
            proteins = proteins.drop(columns="gene")

    peptides = (
        long.groupby(["accession", "peptide"], sort=False)
        .agg(PSMs=("file", "size"), files=("file", "nunique"),
             best_q=("q_value", "min"), charges=("charge", lambda c: ",".join(map(str, sorted(set(c))))),
             shared=("shared", "first"))
        .reset_index()
        .rename(columns={"best_q": "best q-value", "files": "found in files"})
        .sort_values(["accession", "PSMs"], ascending=[True, False], ignore_index=True)
    )
    psms_long = long[["accession", "peptide", "file", "scan", "rt", "mz", "charge", "score", "q_value"]].rename(
        columns={"rt": "RT (s)", "mz": "m/z", "q_value": "q-value"}
    ).reset_index(drop=True)
    return proteins, peptides, psms_long, coverage


def coverage_html(sequence: str, mask: np.ndarray, block: int = 10, per_row: int = 50) -> str:
    """The protein sequence in blocks of ``block`` residues, covered residues highlighted."""
    lines = []
    for row_start in range(0, len(sequence), per_row):
        blocks = []
        for start in range(row_start, min(row_start + per_row, len(sequence)), block):
            chars = []
            for i in range(start, min(start + block, len(sequence))):
                residue = html.escape(sequence[i])
                chars.append(f'<span class="cov">{residue}</span>' if mask[i] else residue)
            blocks.append("".join(chars))
        lines.append(f'<span class="pos">{row_start + 1:>5}</span> ' + " ".join(blocks))
    return (
        "<style>.covmap{font-family:monospace;font-size:14px;line-height:1.7;white-space:pre}"
        ".covmap .cov{background:#d7f5ea;color:#065f46;font-weight:600}"
        ".covmap .pos{color:#a1a1aa}</style>"
        '<div class="covmap">' + "\n".join(lines) + "</div>"
    )

