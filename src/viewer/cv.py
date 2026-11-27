"""Coefficient of variation (CV) of protein abundances within each group.

Streamlit-free. The CV of a protein in a group is the standard deviation of its
measured abundances in the group's samples divided by their mean, in percent,
on the linear scale. Replicates measured precisely give small CVs; the
CVs of all proteins of a group show how reproducible it is.
"""

import pandas as pd

MIN_VALUES = 2  # a CV needs at least two measured values
LOW_CV = 20.0  # percent; proteins below it count as precisely measured


def group_cvs(mat: pd.DataFrame, sample_group_map: dict, min_values: int = MIN_VALUES) -> pd.DataFrame:
    """``protein``, ``group``, ``cv`` (percent) for every protein measured in at
    least ``min_values`` samples of a group.

    ``mat`` is protein x sample on the linear scale; zeros and missing values
    are left out. Groups with fewer than ``min_values`` samples give nothing.
    """
    values = mat.apply(pd.to_numeric, errors="coerce")
    values = values.where(values > 0)
    frames = []
    for group in sorted({g for s, g in sample_group_map.items() if s in values.columns}):
        cols = [s for s in values.columns if sample_group_map.get(s) == group]
        if len(cols) < min_values:
            continue
        sub = values[cols]
        counted = sub.notna().sum(axis=1) >= min_values
        sub = sub[counted]
        cv = sub.std(axis=1, ddof=1) / sub.mean(axis=1) * 100
        frames.append(pd.DataFrame({"protein": sub.index, "group": group, "cv": cv.to_numpy()}))
    if not frames:
        return pd.DataFrame({"protein": pd.Series(dtype=object), "group": pd.Series(dtype=object),
                             "cv": pd.Series(dtype=float)})
    return pd.concat(frames, ignore_index=True)


def cv_summary(cvs: pd.DataFrame, sample_group_map: dict, low: float = LOW_CV) -> pd.DataFrame:
    """Per group: ``samples``, ``proteins`` with a CV, ``median`` CV and the
    percentage of proteins below ``low``."""
    rows = []
    for group, sub in cvs.groupby("group", sort=True):
        rows.append({
            "group": group,
            "samples": sum(1 for g in sample_group_map.values() if g == group),
            "proteins": len(sub),
            "median": float(sub["cv"].median()),
            "below": float((sub["cv"] < low).mean() * 100),
        })
    return pd.DataFrame(rows, columns=["group", "samples", "proteins", "median", "below"])


def plot_table(cvs: pd.DataFrame, accession=str) -> pd.DataFrame:
    """The library's contract for the boxes: ``condition``, ``cv`` and the
    protein's ``accession`` (``accession(protein)``), by which a picked protein
    is marked."""
    return pd.DataFrame({
        "condition": cvs["group"].astype(str),
        "cv": cvs["cv"].astype(float),
        "accession": cvs["protein"].map(accession).astype(str),
    })
