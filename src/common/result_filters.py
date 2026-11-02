"""Filters of the PSM and protein tables, without Streamlit.

The tables filter themselves (openms-insight's search box and filter dialog);
this module adds the columns they filter on and the thresholds they offer.
Every table on those pages holds only what passed the run's PSM FDR, so a
q-value threshold can only tighten it: the values offered never exceed the
run's level.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.common.contaminants import is_contaminant
from src.common.fdr_display import percent_label

# Tighter q-value levels offered below the run's PSM FDR, in percent
STANDARD_LEVELS = (5.0, 1.0, 0.5, 0.1, 0.01)


def q_choices(run_fdr_percent: float | None) -> list[float]:
    """q-value levels (percent) a viewer may pick: the run's own first, then tighter ones.

    The run's level is None when its PSM FDR filter was off; then every standard
    level is offered, below 100 % (no threshold).
    """
    top = 100.0 if run_fdr_percent is None else float(run_fdr_percent)
    return [top] + [level for level in STANDARD_LEVELS if level < top - 1e-12]


def q_suggestions(run_fdr_percent: float | None) -> list[dict]:
    """q-value thresholds for the PSM table's filter dialog, the run's own level first.

    Each is ``{"value": q, "label": ...}``, as openms-insight's ``numeric_filters``
    takes them; without a run level (PSM FDR off) the standard levels are offered.
    """
    choices = q_choices(run_fdr_percent)
    suggestions = []
    for i, level in enumerate(choices):
        if i == 0 and run_fdr_percent is None:
            continue  # "no threshold" is the dialog's own "Any"
        label = f"{percent_label(level)} (run)" if i == 0 else percent_label(level)
        suggestions.append({"value": level / 100.0, "label": label})
    return suggestions


def protein_numeric_filters(max_replicates: int) -> dict:
    """``numeric_filters`` of the protein table: min. peptides and, with conditions, replicates."""
    filters = {"peptides": {"op": ">=", "suggestions": [2, 3, 5]}}
    if max_replicates >= 1:
        filters["replicates"] = {"op": ">=", "suggestions": list(range(1, max_replicates + 1))}
    return filters


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
