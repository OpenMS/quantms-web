"""Filters of the filter bars on the PSM and protein pages, without Streamlit.

Every table on those pages holds only what passed the run's PSM FDR, so a
q-value threshold here can only tighten it; the choices offered never exceed
the run's level.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from src.common.contaminants import is_contaminant

# Tighter q-value levels offered below the run's PSM FDR, in percent
STANDARD_LEVELS = (5.0, 1.0, 0.5, 0.1, 0.01)


def q_choices(run_fdr_percent: float | None) -> list[float]:
    """q-value levels (percent) a viewer may pick: the run's own first, then tighter ones.

    The run's level is None when its PSM FDR filter was off; then every standard
    level is offered, below 100 % (no threshold).
    """
    top = 100.0 if run_fdr_percent is None else float(run_fdr_percent)
    return [top] + [level for level in STANDARD_LEVELS if level < top - 1e-12]


def text_mask(df: pd.DataFrame, columns, query: str) -> pd.Series:
    """Rows where any of ``columns`` contains ``query`` (case-insensitive, literal).

    Several words must all match, each in any of the columns.
    """
    mask = pd.Series(True, index=df.index)
    for word in (query or "").split():
        pattern = re.escape(word)
        hit = pd.Series(False, index=df.index)
        for column in columns:
            if column in df.columns:
                hit |= df[column].astype(str).str.contains(pattern, case=False, regex=True, na=False)
        mask &= hit
    return mask


def filter_psms_by_q(psms: pd.DataFrame, q_max: float | None) -> pd.DataFrame:
    """PSMs with q-value <= q_max (a fraction); PSMs without a q-value are kept."""
    if q_max is None:
        return psms
    q = psms["q_value"]
    return psms[q.isna() | (q <= q_max + 1e-12)]


def file_conditions(files, group_of) -> dict:
    """File -> condition, for the files that have one (``group_of(file)`` gives it or "")."""
    return {f: c for f in files if (c := str(group_of(f) or "").strip())}


def replicates_found(psms: pd.DataFrame, conditions: dict, sample_of: dict | None = None) -> pd.Series:
    """Per accession, the most files of one condition with a PSM of it.

    ``psms`` is long format (one row per PSM and accession, column ``file``);
    files without a condition are left out. 0 for a protein seen only there.
    With ``sample_of`` (file -> sample, for fractionated runs) samples are
    counted instead of files, so the fractions of one sample count once.
    """
    if psms.empty or not conditions:
        return pd.Series(dtype="int64")
    seen = psms[["accession", "file"]].drop_duplicates()
    seen = seen.assign(
        condition=seen["file"].map(conditions),
        unit=seen["file"].map(lambda f: (sample_of or {}).get(f, f)),
    ).dropna(subset=["condition"])
    if seen.empty:
        return pd.Series(dtype="int64")
    return seen.groupby(["accession", "condition"])["unit"].nunique().groupby(level=0).max()


def add_protein_flags(
    proteins: pd.DataFrame, psms_long: pd.DataFrame, conditions: dict, sample_of: dict | None = None,
) -> pd.DataFrame:
    """Add ``contaminant`` and, with conditions set, ``replicates`` (see replicates_found)."""
    proteins = proteins.copy()
    proteins["contaminant"] = proteins["accession"].map(is_contaminant).astype(bool)
    if conditions:
        counts = replicates_found(psms_long, conditions, sample_of)
        replicates = proteins["accession"].map(counts).fillna(0).astype(int)
        # next to "found in", ahead of the wide per-file columns
        at = proteins.columns.get_loc("found in") + 1 if "found in" in proteins.columns else len(proteins.columns)
        proteins.insert(at, "replicates", replicates)
    return proteins


def filter_proteins(
    proteins: pd.DataFrame,
    query: str = "",
    min_peptides: int = 1,
    min_replicates: int = 0,
) -> pd.DataFrame:
    """Proteins matching the search, with at least ``min_peptides`` peptides and, if
    the table has a ``replicates`` column, found in ``min_replicates`` files of one condition."""
    keep = text_mask(proteins, ("accession", "gene", "description"), query)
    keep &= proteins["peptides"] >= int(min_peptides)
    if min_replicates and "replicates" in proteins.columns:
        keep &= proteins["replicates"] >= int(min_replicates)
    return proteins[keep]


def max_replicates(conditions: dict, sample_of: dict | None = None) -> int:
    """Files (or samples, with ``sample_of``) in the largest condition: the highest
    "found in" threshold that means anything."""
    if not conditions:
        return 0
    units = pd.DataFrame({
        "condition": list(conditions.values()),
        "unit": [(sample_of or {}).get(f, f) for f in conditions],
    })
    return int(units.groupby("condition")["unit"].nunique().max())


def q_value_per_hit(peptides) -> np.ndarray:
    """q-value (MS:1001491, else "q-value") of every hit, in parse_idxml's id_idx order."""
    values = []
    for pep in peptides:
        for hit in pep.getHits():
            q = np.nan
            for key in ("MS:1001491", "q-value"):
                if hit.metaValueExists(key):
                    q = float(hit.getMetaValue(key))
                    break
            values.append(q)
    return np.asarray(values, dtype="float64")
