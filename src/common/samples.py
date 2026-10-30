"""Per-sample counts of a finished run: PSMs, peptides and proteins, by condition.

Reads only what the workflow already wrote: the IDFilter idXMLs (identifications
at the run's PSM FDR) and the quantification result. It makes no Streamlit
calls, so the page caches it and tests call it.
"""
from pathlib import Path

import pandas as pd

from src.common.results_helpers import extract_filename_from_idxml, load_idxml

UNASSIGNED = "Unassigned"

SAMPLE_COLUMNS = [
    "Sample",
    "Condition",
    "Target PSMs",
    "Peptides",
    "Proteins",
    "Quantified proteins",
    "Quantified peptides",
]

# Columns the bar chart can show, in the order the page offers them
COUNT_COLUMNS = SAMPLE_COLUMNS[2:]


def _is_decoy(hit) -> bool:
    if not hit.metaValueExists("target_decoy"):
        return False
    value = hit.getMetaValue("target_decoy")
    if isinstance(value, bytes):
        value = value.decode()
    return str(value).startswith("decoy")


def identification_counts(idxml_path) -> dict:
    """Target PSMs, distinct peptides and distinct proteins of one idXML.

    Only the top hit of each spectrum counts, and decoy hits are left out.
    Proteins are the accessions the target peptides map to.
    """
    _, peptides = load_idxml(idxml_path)
    psms = 0
    sequences, proteins = set(), set()
    for pep in peptides:
        hits = pep.getHits()
        if not hits or _is_decoy(hits[0]):
            continue
        hit = hits[0]
        psms += 1
        sequences.add(hit.getSequence().toString())
        proteins.update(ev.getProteinAccession() for ev in hit.getPeptideEvidences())
    return {"Target PSMs": psms, "Peptides": len(sequences), "Proteins": len(proteins)}


def lfq_quantified_counts(quant_csv) -> dict:
    """Quantified proteins and peptides per input file name, from ProteomicsLFQ's CSV.

    A protein (group) or peptide counts for a file when its intensity there is above 0.
    """
    try:
        df = pd.read_csv(quant_csv, usecols=["ProteinName", "PeptideSequence", "Reference", "Intensity"])
    except (OSError, ValueError):
        return {}
    df = df[pd.to_numeric(df["Intensity"], errors="coerce") > 0]
    counts = df.groupby("Reference").agg(
        proteins=("ProteinName", "nunique"), peptides=("PeptideSequence", "nunique")
    )
    return {
        Path(str(ref)).name: {"Quantified proteins": int(row.proteins), "Quantified peptides": int(row.peptides)}
        for ref, row in counts.iterrows()
    }


def _condition(group_map: dict, *names) -> str:
    for name in names:
        value = group_map.get(name)
        if value is not None and str(value).strip():
            return str(value).strip()
    return UNASSIGNED


def _order(df: pd.DataFrame) -> pd.DataFrame:
    """Samples grouped by condition (unassigned last), in their original order within one."""
    conditions = sorted(c for c in df["Condition"].unique() if c != UNASSIGNED)
    if UNASSIGNED in set(df["Condition"]):
        conditions.append(UNASSIGNED)
    df = df.copy()
    df["Condition"] = pd.Categorical(df["Condition"], categories=conditions, ordered=True)
    df = df.sort_values("Condition", kind="stable").reset_index(drop=True)
    df["Condition"] = df["Condition"].astype(str)
    return df


def lfq_samples(workflow_dir, mzml_files, group_map: dict) -> pd.DataFrame:
    """One row per input file (SAMPLE_COLUMNS); files the run did not reach show empty cells.

    ``group_map`` maps an input file name to its condition, as set on the Configure page.
    """
    results = Path(workflow_dir) / "results"
    names = [Path(f).name for f in mzml_files]
    filter_dir = results / "psm_filter"
    idxmls = {}
    if filter_dir.is_dir():
        idxmls = {extract_filename_from_idxml(p, names): p for p in sorted(filter_dir.glob("*.idXML"))}
    quant_csvs = sorted((results / "quant_results").glob("*.csv")) if (results / "quant_results").is_dir() else []
    quantified = lfq_quantified_counts(quant_csvs[0]) if quant_csvs else {}

    rows = []
    for name in names:
        ids = identification_counts(idxmls[name]) if name in idxmls else {}
        rows.append({
            "Sample": Path(name).stem,
            "Condition": _condition(group_map, name, Path(name).stem),
            **{col: ids.get(col) for col in ("Target PSMs", "Peptides", "Proteins")},
            **{col: quantified.get(name, {}).get(col) for col in ("Quantified proteins", "Quantified peptides")},
        })
    return _order(pd.DataFrame(rows, columns=SAMPLE_COLUMNS))


def channel_samples(pivot_df: pd.DataFrame, sample_columns, sample_group_map: dict) -> pd.DataFrame:
    """One row per isobaric channel: the proteins quantified in it (intensity above 0).

    Identifications are per file, not per channel, so those columns stay empty.
    """
    rows = []
    for col in sample_columns:
        values = pd.to_numeric(pivot_df[col], errors="coerce")
        rows.append({
            "Sample": col,
            "Condition": _condition(sample_group_map, col),
            "Quantified proteins": int((values > 0).sum()),
        })
    return _order(pd.DataFrame(rows, columns=SAMPLE_COLUMNS))


def condition_summary(samples: pd.DataFrame) -> pd.DataFrame:
    """Per condition: the number of samples and the mean of each count over them."""
    counts = [c for c in COUNT_COLUMNS if samples[c].notna().any()]
    numeric = samples[["Condition"] + counts].copy()
    numeric[counts] = numeric[counts].apply(pd.to_numeric, errors="coerce")
    summary = numeric.groupby("Condition", sort=False).agg(
        Samples=("Condition", "size"), **{f"Mean {c}": (c, "mean") for c in counts}
    )
    return summary.round(1).reset_index()
