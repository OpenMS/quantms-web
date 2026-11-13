"""Diagnostic plots that help choose a setting on each step of the group comparison.

``filter_threshold_curve`` returns a Plotly figure. The others return tidy
tables for an openms-insight ``SeriesPlot`` (``show_series``): a position
column, a value column and, where several series are drawn, a series column.
"""
from typing import Callable

import numpy as np
import pandas as pd
import plotly.graph_objects as go


def _log2(values: pd.DataFrame) -> pd.DataFrame:
    """log2 intensities with zeros and non-positive values treated as missing."""
    values = values.apply(pd.to_numeric, errors="coerce")
    return np.log2(values.where(values > 0))


def filter_threshold_curve(
    run_filter: Callable[[float], int],
    thresholds: list[float],
    current: float,
    total: int,
    x_label: str,
) -> go.Figure:
    """Proteins kept by a filter across its threshold range.

    ``run_filter(threshold)`` returns the number of proteins the filter keeps.
    A plateau means the threshold is not sensitive there; a steep drop shows
    where the filter starts removing many proteins.
    """
    kept = [run_filter(t) for t in thresholds]
    fig = go.Figure(go.Scatter(x=thresholds, y=kept, mode="lines+markers", name="Proteins kept"))
    fig.add_hline(y=total, line_dash="dot", line_color="grey", annotation_text=f"All proteins ({total})")
    fig.add_vline(x=current, line_dash="dash", line_color="#E74C3C", annotation_text="Current")
    fig.update_layout(
        title="Proteins kept at each threshold",
        xaxis_title=x_label,
        yaxis_title="Proteins kept",
        yaxis_rangemode="tozero",
        showlegend=False,
    )
    return fig


def _bin_centres(edges: np.ndarray) -> np.ndarray:
    return np.round((edges[:-1] + edges[1:]) / 2, 4)


def _histogram_lines(values: dict[str, pd.Series], bins: int, density: bool) -> pd.DataFrame:
    """Per series, the counts (or density) of ``values`` over shared bins, as x, y, series."""
    finite = pd.concat([v for v in values.values()]).replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return pd.DataFrame({"x": [], "y": [], "series": []})
    edges = np.histogram_bin_edges(finite.to_numpy(dtype=float), bins=bins)
    frames = []
    for name, series in values.items():
        data = series.replace([np.inf, -np.inf], np.nan).dropna().to_numpy(dtype=float)
        if not len(data):
            continue
        counts, _ = np.histogram(data, bins=edges, density=density)
        frames.append(pd.DataFrame({"x": _bin_centres(edges), "y": counts, "series": name}))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame({"x": [], "y": [], "series": []})


MISSING_SERIES = {"complete": "Complete", "missing": "With missing values"}
MISSING_COLORS = {"complete": "#3498DB", "missing": "#E67E22"}


def missingness_vs_intensity(df: pd.DataFrame, sample_cols: list[str]) -> tuple[pd.DataFrame, int]:
    """Density of the mean log2 intensity of proteins with and without missing values.

    If proteins with missing values sit clearly lower than complete ones,
    values are missing because they fall below the detection limit (MNAR);
    overlapping distributions point to random dropout (MAR). Returns the
    table (x, y, series) and the number of proteins with a gap.
    """
    log_values = _log2(df[sample_cols])
    has_missing = log_values.isna().any(axis=1)
    means = log_values.mean(axis=1)
    table = _histogram_lines(
        {"complete": means[~has_missing], "missing": means[has_missing]}, bins=40, density=True
    )
    return table, int(has_missing.sum())


IMPUTED_SERIES = {"observed": "Observed", "imputed": "Filled in"}
IMPUTED_COLORS = {"observed": "#3498DB", "imputed": "#E74C3C"}


def imputed_value_preview(before: pd.DataFrame, after: pd.DataFrame, sample_cols: list[str]) -> pd.DataFrame:
    """Counts of observed log2 values and of the values imputation fills in, as x, y, series."""
    before_log = _log2(before[sample_cols])
    after_log = _log2(after[sample_cols])
    observed = before_log.stack()
    imputed = after_log.where(before_log.isna()).stack()
    return _histogram_lines({"observed": observed, "imputed": imputed}, bins=60, density=False)


STAGE_SERIES = {"before": "Before", "after": "After"}
STAGE_COLORS = {"before": "#95A5A6", "after": "#3498DB"}


def sample_distributions(
    before: pd.DataFrame,
    after: pd.DataFrame,
    sample_cols: list[str],
    after_is_log: bool,
    points: int = 201,
) -> pd.DataFrame:
    """Per sample, log2 intensities before and after normalization, as Sample, Intensity, Stage.

    Well-normalized samples have aligned medians; a sample that stays offset
    after normalization is worth checking. Untransformed tables are shown on
    a log2 scale. Each sample and stage keeps ``points`` evenly spaced
    quantiles of its values, which draw the same box as all of them.
    """
    frames = []
    quantiles = np.linspace(0, 1, points)
    for stage, table, is_log in [("before", before, False), ("after", after, after_is_log)]:
        values = table[sample_cols].apply(pd.to_numeric, errors="coerce")
        if not is_log:
            values = _log2(values)
        for sample in sample_cols:
            column = values[sample].replace([np.inf, -np.inf], np.nan).dropna()
            if column.empty:
                continue
            kept = column.to_numpy(dtype=float) if len(column) <= points else column.quantile(quantiles).to_numpy()
            frames.append(pd.DataFrame({"Sample": sample, "Intensity": kept, "Stage": stage}))
    if not frames:
        return pd.DataFrame({"Sample": [], "Intensity": [], "Stage": []})
    return pd.concat(frames, ignore_index=True)


def pvalue_histogram(statistics_df: pd.DataFrame, bins: int = 20) -> tuple[pd.DataFrame, float]:
    """Proteins per p-value bin (x = bin centre, y = proteins) and the count expected per bin.

    A flat histogram with a peak near 0 is what a sound test looks like. A
    U-shape or a peak near 1 suggests the input scale or the chosen test does
    not fit the data. The expected count is what each bin holds when no
    protein changes.
    """
    pvalues = pd.to_numeric(statistics_df["p-value"], errors="coerce").dropna()
    edges = np.linspace(0, 1, bins + 1)
    counts, _ = np.histogram(pvalues.clip(0, 1).to_numpy(dtype=float), bins=edges)
    table = pd.DataFrame({"p-value": _bin_centres(edges), "Proteins": counts})
    return table, len(pvalues) / bins
