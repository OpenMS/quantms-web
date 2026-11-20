"""Coefficient of variation (CV) of protein abundances within each group.

Streamlit-free. The CV of a protein in a group is the standard deviation of its
measured abundances in the group's samples divided by their mean, in percent,
on the linear scale. Replicates measured precisely give small CVs; the
distribution over all proteins shows how reproducible each group is.
"""

import numpy as np
import pandas as pd

MIN_VALUES = 2  # a CV needs at least two measured values
LOW_CV = 20.0  # percent; proteins below it count as precisely measured
BIN_WIDTH = 4.0  # percent


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


def cv_upper(cvs: pd.DataFrame) -> float:
    """The x range of the distribution: the 99th percentile, rounded up to 10, 50 to 200 %."""
    if cvs.empty:
        return 100.0
    top = float(np.nanpercentile(cvs["cv"], 99))
    return float(min(200.0, max(50.0, np.ceil(top / 10) * 10)))


def cv_distribution(cvs: pd.DataFrame, bin_width: float = BIN_WIDTH, upper: float | None = None) -> pd.DataFrame:
    """Tidy ``cv``, ``proteins`` (percent of the group's proteins), ``group``:
    one point per bin centre, for a line per group. CVs above ``upper`` go in
    the last bin, so every group's curve sums to 100 %."""
    upper = cv_upper(cvs) if upper is None else upper
    edges = np.arange(0.0, upper + bin_width, bin_width)
    centres = (edges[:-1] + edges[1:]) / 2
    rows = []
    for group, sub in cvs.groupby("group", sort=True):
        counts, _ = np.histogram(np.clip(sub["cv"].to_numpy(), 0, upper - bin_width / 2), bins=edges)
        share = counts / max(len(sub), 1) * 100
        rows.append(pd.DataFrame({"cv": centres, "proteins": share, "group": group}))
    if not rows:
        return pd.DataFrame({"cv": [], "proteins": [], "group": []})
    return pd.concat(rows, ignore_index=True)


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
