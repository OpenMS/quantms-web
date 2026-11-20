"""Coefficient of variation within groups (src/common/cv.py)."""

import numpy as np
import pandas as pd
import pytest

from src.common.cv import cv_distribution, cv_summary, cv_upper, group_cvs

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


def test_distribution_sums_to_100_percent_per_group():
    cvs = pd.DataFrame({"protein": list("pqrstu"), "group": ["A"] * 3 + ["B"] * 3,
                        "cv": [1.0, 15.0, 500.0, 3.0, 3.0, 41.0]})
    curves = cv_distribution(cvs, bin_width=2, upper=50)
    for _, sub in curves.groupby("group"):
        assert sub["proteins"].sum() == pytest.approx(100)
    # 500 % goes in the last bin
    a = curves[curves["group"] == "A"].set_index("cv")["proteins"]
    assert a.iloc[-1] == pytest.approx(100 / 3)
    assert curves["cv"].max() < 50


def test_upper_end_of_the_axis():
    assert cv_upper(pd.DataFrame({"cv": [5.0] * 100})) == 50
    assert cv_upper(pd.DataFrame({"cv": np.linspace(0, 123, 1000)})) == 130
    assert cv_upper(pd.DataFrame({"cv": [1000.0] * 10})) == 200


def test_summary_per_group():
    summary = cv_summary(group_cvs(_matrix(), GROUPS), GROUPS).set_index("group")
    assert summary.loc["A", "samples"] == 3
    assert summary.loc["A", "proteins"] == 3
    assert summary.loc["B", "proteins"] == 4
    assert summary.loc["B", "below"] == pytest.approx(50)
