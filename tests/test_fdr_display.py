"""FDR wording shared by Configure and the results pages."""

from src.common.fdr_display import (
    FDR_QUICK_CHOICES,
    fdr_applied,
    fdr_label,
    matching_quick_choice,
    percent_label,
    q_label,
)


def test_percent_and_q_labels_agree():
    assert percent_label(1.0) == "1 %"
    assert q_label(1.0) == "q ≤ 0.01"
    assert percent_label(0.1) == "0.1 %"
    assert q_label(0.1) == "q ≤ 0.001"
    assert q_label(5) == "q ≤ 0.05"
    # no float noise from the division
    assert q_label(1.5) == "q ≤ 0.015"
    assert q_label(0.3) == "q ≤ 0.003"


def test_100_percent_is_not_applied():
    assert not fdr_applied(100.0)
    assert not fdr_applied(None)
    assert fdr_applied(99.9)
    assert fdr_label(100) == "not applied"
    assert fdr_label(None) == "not applied"
    assert fdr_label(1) == "1 % (q ≤ 0.01)"


def test_quick_choices_match_typed_values():
    assert FDR_QUICK_CHOICES == (0.1, 1.0, 5.0)
    assert matching_quick_choice(1) == 1.0
    assert matching_quick_choice(0.1) == 0.1
    assert matching_quick_choice(2.0) is None
    assert matching_quick_choice(None) is None
