"""How an FDR level reads on screen: "1 %" next to "q ≤ 0.01", or "not applied".

Streamlit-free, so the workflow, the results pages and the tests share one wording.
The PSM FDR is stored as a percentage ("psm-fdr-percent"); 100 % switches the filter off.
"""

from __future__ import annotations

# Quick choices offered next to the number field, in percent.
FDR_QUICK_CHOICES = (0.1, 1.0, 5.0)


def _num(value: float) -> str:
    """A number without float noise or trailing zeros (0.1, 1, 0.015)."""
    return f"{round(float(value), 10):g}"


def fdr_applied(percent: float | None) -> bool:
    """True when an FDR filter at this level removes anything (below 100 %)."""
    return percent is not None and float(percent) < 100.0


def percent_label(percent: float) -> str:
    """'1 %' for 1.0."""
    return f"{_num(percent)} %"


def q_label(percent: float) -> str:
    """'q ≤ 0.01' for 1 %."""
    return f"q ≤ {_num(float(percent) / 100.0)}"


def fdr_label(percent: float | None) -> str:
    """'1 % (q ≤ 0.01)', or 'not applied' when the filter was off or skipped."""
    if not fdr_applied(percent):
        return "not applied"
    return f"{percent_label(percent)} (q ≤ {_num(float(percent) / 100.0)})"


def matching_quick_choice(percent: float | None) -> float | None:
    """The quick choice equal to this level, so the chip shows as selected."""
    if percent is None:
        return None
    for choice in FDR_QUICK_CHOICES:
        if abs(float(percent) - choice) < 1e-9:
            return choice
    return None
