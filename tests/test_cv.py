"""Coefficient of variation within groups (src/common/cv.py)."""

import numpy as np
import pandas as pd
import pytest

from src.viewer.cv import cv_summary, group_cvs, plot_table

GROUPS = {"a1": "A", "a2": "A", "a3": "A", "b1": "B", "b2": "B", "c1": "C"}


def _matrix():
    return pd.DataFrame(
        {
            "a1": [100.0, 10.0, 0.0, 5.0],
            "a2": [100.0, 20.0, 7.0, np.nan],
            "a3": [100.0, 30.0, 9.0, np.nan],
            "b1": [50.0, 1.0, 1.0, 2.0],
            "b2": [150.0, 1.0, 3.0, 2.0],
            "c1": [1.0, 2.0, 3.0, 4.0],
        },
        index=["P1", "P2", "P3", "P4"],
    )


def test_cv_per_group_uses_measured_values_only():
    cvs = group_cvs(_matrix(), GROUPS).set_index(["group", "protein"])["cv"]
    assert cvs["A", "P1"] == pytest.approx(0.0)
    assert cvs["A", "P2"] == pytest.approx(np.std([10, 20, 30], ddof=1) / 20 * 100)
    # A zero is not measured: P3 counts its two other values in A
    assert cvs["A", "P3"] == pytest.approx(np.std([7, 9], ddof=1) / 8 * 100)
    # One measured value in A gives no CV
    assert ("A", "P4") not in cvs.index
    assert cvs["B", "P1"] == pytest.approx(np.std([50, 150], ddof=1) / 100 * 100)
    # A group of one sample has no CV at all
    assert "C" not in cvs.index.get_level_values("group")


def test_no_group_with_two_samples_gives_nothing():
    assert group_cvs(_matrix(), {"a1": "A", "b1": "B"}).empty


def test_plot_table_names_the_lead_accession():
    cvs = pd.DataFrame({"protein": ["P1;P2", "Q1"], "group": ["A", "B"], "cv": [12.5, 30.0]})
    table = plot_table(cvs, lambda p: p.split(";")[0])
    assert list(table.columns) == ["condition", "cv", "accession"]
    assert table["accession"].tolist() == ["P1", "Q1"]
    assert table["condition"].tolist() == ["A", "B"]


def test_summary_per_group():
    summary = cv_summary(group_cvs(_matrix(), GROUPS), GROUPS).set_index("group")
    assert summary.loc["A", "samples"] == 3
    assert summary.loc["A", "proteins"] == 3
    assert summary.loc["B", "proteins"] == 4
    assert summary.loc["B", "below"] == pytest.approx(50)
