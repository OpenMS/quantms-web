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

# run and psm_id locate the PSM in the run's spectrum viewer cache (see parse_idxml)
PSM_COLUMNS = [
    "file", "scan", "rt", "mz", "charge", "peptide", "sequence", "score", "q_value", "proteins", "run", "psm_id",
]


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
        # parse_idxml numbers every hit of the run in this order; the top hit's
        # number is its id in the run's SequenceView cache
        n_hits = 0
        for pep in load_idxml(path)[1]:
            hits = pep.getHits()
            psm_id = n_hits
            n_hits += len(hits)
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
                "run": Path(path).stem,
                "psm_id": psm_id,
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


def build_tables(
    psms: pd.DataFrame,
    fasta: dict[str, tuple[str, str]],
    files: list[str],
    sample_of_file: dict[str, str] | None = None,
):
    """Proteins, peptide-to-protein and PSM-to-protein tables.

    With ``sample_of_file`` (file name to sample name, for a run whose samples
    were split into fractions), the per-file counts are per sample instead.

    Peptides and PSMs are in long format (one row per protein they map to), so
    selecting a protein filters them by its accession.

    Returns:
        (proteins, peptides, psms_long, coverage) where ``coverage`` maps an
        accession to its residue mask.
    """
    long = psms.explode("proteins").rename(columns={"proteins": "accession"})
    long["shared"] = psms["proteins"].str.len().reindex(long.index).gt(1)
    if sample_of_file:
        long["sample"] = long["file"].map(lambda f: sample_of_file.get(f, Path(f).stem))
        columns = list(dict.fromkeys(sample_of_file.get(f, Path(f).stem) for f in files))
    else:
        long["sample"] = long["file"]
        columns = list(files)

    rows, coverage = [], {}
    for accession, group in long.groupby("accession", sort=False):
        description, sequence = fasta.get(accession, ("", ""))
        mask = coverage_mask(sequence, group["sequence"])
        coverage[accession] = mask
        per_file = group.groupby("sample").size()
        row = {
            "accession": accession,
            "gene": gene_name(description),
            "description": protein_name(description),
            "coverage (%)": round(100 * mask.mean(), 1) if len(mask) else np.nan,
            "peptides": group["sequence"].nunique(),
            "unique peptides": group.loc[~group["shared"], "sequence"].nunique(),
            "PSMs": len(group),
            "found in": f"{per_file.size} of {len(columns)}",
            "length": len(sequence) or np.nan,
        }
        for column in columns:
            label = column if sample_of_file else Path(column).stem
            row[f"PSMs {label}"] = int(per_file.get(column, 0))
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
        .agg(PSMs=("file", "size"), files=("sample", "nunique"),
             best_q=("q_value", "min"), charges=("charge", lambda c: ",".join(map(str, sorted(set(c))))),
             shared=("shared", "first"))
        .reset_index()
        .rename(columns={"best_q": "best q-value", "files": "found in samples" if sample_of_file else "found in files"})
        .sort_values(["accession", "PSMs"], ascending=[True, False], ignore_index=True)
    )
    psms_long = long[["accession", "peptide", "file", "scan", "rt", "mz", "charge", "score", "q_value"]].rename(
        columns={"rt": "RT (s)", "mz": "m/z", "q_value": "q-value"}
    ).reset_index(drop=True)
    if {"run", "psm_id"} <= set(long.columns):
        # One value per PSM, for selecting it: "<run>:<psm_id>"
        psms_long["psm"] = (long["run"].astype(str) + ":" + long["psm_id"].astype(str)).to_numpy()
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



def _peptide_modifications(peptide: str, offset: int) -> list[dict]:
    """Modified residues of ``peptide`` (an OpenMS string), as positions in the
    protein (``offset``: the 1-based start of the peptide in it) and names."""
    try:
        from pyopenms import AASequence

        seq = AASequence.fromString(peptide)
    except Exception:  # pylint: disable=broad-except
        return []
    mods = []
    if seq.hasNTerminalModification():
        mods.append({"position": offset, "token": seq.getNTerminalModificationName()})
    for i in range(seq.size()):
        residue = seq.getResidue(i)
        if residue.isModified():
            mods.append({"position": offset + i, "token": residue.getModificationName()})
    if seq.hasCTerminalModification():
        mods.append({"position": offset + seq.size() - 1, "token": seq.getCTerminalModificationName()})
    return mods


COVERAGE_COLUMNS = ("match_id", "reference", "track", "start", "end", "sequence", "count",
                    "mod_positions", "mod_labels")


def coverage_rows(psms_long: pd.DataFrame, fasta: dict, track_of_file: dict | None = None,
                  default_track: str = "All files"):
    """Where each peptide lies on its protein, one row per protein, peptide and
    track (the condition of the PSM's file, else ``default_track``), counted
    in PSMs: the column contract of the sequence coverage panel
    (``src/viz/identification/sequence_coverage.py``).

    ``match_id`` is the peptide (as the peptide table names it), so picking a
    bar picks the peptide; ``reference`` is the protein accession. A peptide
    found more than once in its protein is drawn at its first place; one not
    found in the FASTA sequence (or a protein without one) is left out.
    """
    import polars as pl

    schema = {"match_id": pl.Utf8, "reference": pl.Utf8, "track": pl.Utf8, "start": pl.Int64,
              "end": pl.Int64, "sequence": pl.Utf8, "count": pl.Int64,
              "mod_positions": pl.List(pl.Int64), "mod_labels": pl.List(pl.Utf8)}
    if psms_long.empty:
        return pl.DataFrame(schema=schema)
    tracks = psms_long["file"].map(lambda f: (track_of_file or {}).get(f) or default_track)
    counts = psms_long.assign(track=tracks).groupby(["accession", "peptide", "track"], sort=False).size()
    records, places = [], {}
    for (accession, peptide, track), count in counts.items():
        key = (accession, peptide)
        if key not in places:
            places[key] = None
            sequence = fasta.get(accession, ("", ""))[1]
            plain = _unmodified(peptide)
            at = sequence.find(plain) if sequence and plain else -1
            if at >= 0:
                places[key] = (at + 1, at + len(plain), plain, _peptide_modifications(peptide, at + 1))
        place = places[key]
        if place is None:
            continue
        start, end, plain, mods = place
        records.append({
            "match_id": peptide, "reference": accession, "track": track, "start": start, "end": end,
            "sequence": plain, "count": int(count),
            "mod_positions": [m["position"] for m in mods], "mod_labels": [m["token"] for m in mods],
        })
    return pl.DataFrame(records, schema=schema)


def _unmodified(peptide: str) -> str:
    """``peptide`` without its modifications (``.(Acetyl)PEPM(Oxidation)K`` -> ``PEPMK``)."""
    try:
        from pyopenms import AASequence

        return AASequence.fromString(peptide).toUnmodifiedString()
    except Exception:  # pylint: disable=broad-except
        return re.sub(r"\([^)]*\)|\[[^\]]*\]|\.", "", peptide)
