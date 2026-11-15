"""The filter step's threshold curve (a Plotly figure). The plots of the other
steps of the group comparison are the template's (``src/viz/comparison_steps.py``).
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
