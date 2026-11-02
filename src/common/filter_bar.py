"""Filter bars above the PSM and protein tables.

The widgets keep their values per page in session state (not workspace
parameters: they change what is shown, not what the analysis did).
"""
import streamlit as st

from src.common.fdr_display import percent_label, q_label
from src.common.result_filters import q_choices


def _q_select(key: str, run_fdr_percent: float | None, column) -> float | None:
    """q-value threshold; returns a fraction, or None for the run's own level."""
    choices = q_choices(run_fdr_percent)

    def label(level: float) -> str:
        if level == choices[0]:
            if run_fdr_percent is None:
                return "any (no PSM FDR in this run)"
            return f"{percent_label(level)} · {q_label(level)} (run)"
        return f"{percent_label(level)} · {q_label(level)}"

    level = column.selectbox(
        "PSM q-value",
        choices,
        format_func=label,
        key=f"{key}-q",
        help="Show only PSMs that pass a stricter false discovery rate than the run used. "
             "The run's own level is the loosest choice: these results hold nothing above it.",
    )
    return None if level == choices[0] else level / 100.0


def psm_filter_bar(key: str, run_fdr_percent: float | None) -> tuple[str, float | None]:
    """Search and q-value threshold for a PSM table. Returns (query, q_max)."""
    with st.container(border=True):
        c1, c2 = st.columns([3, 2])
        query = c1.text_input(
            "Search",
            key=f"{key}-search",
            placeholder="Peptide sequence or protein accession, e.g. LVNELTEFAK or P02768",
            help="Rows whose sequence (in any modification form) or proteins contain every word. "
                 "Not case sensitive.",
        )
        q_max = _q_select(key, run_fdr_percent, c2)
        st.caption("Hide contaminants with the switch in the table's title row.")
    return query, q_max


def protein_filter_bar(
    key: str, run_fdr_percent: float | None, max_replicates: int
) -> tuple[str, float | None, int, int]:
    """Search, q-value, minimum peptides and replicates for the protein table.

    Returns (query, q_max, min_peptides, min_replicates); max_replicates is the
    size of the largest condition, 0 when no conditions are set.
    """
    with st.container(border=True):
        c1, c2, c3, c4 = st.columns([3, 2, 1.3, 1.7])
        query = c1.text_input(
            "Search",
            key=f"{key}-search",
            placeholder="Accession, gene or protein name",
            help="Proteins whose accession, gene or description contains every word. "
                 "Not case sensitive.",
        )
        q_max = _q_select(key, run_fdr_percent, c2)
        min_peptides = c3.number_input(
            "Min. peptides",
            min_value=1,
            max_value=50,
            step=1,
            key=f"{key}-min-peptides",
            help="Proteins with at least this many distinct peptides (after the q-value filter). "
                 "2 leaves out one-hit wonders.",
        )
        if max_replicates >= 1:
            min_replicates = c4.number_input(
                "Found in ≥ n replicates",
                min_value=0,
                max_value=max_replicates,
                step=1,
                key=f"{key}-min-replicates",
                help="Proteins with PSMs in at least this many files of one condition "
                     "(any condition). 0 shows every protein.",
            )
        else:
            min_replicates = 0
            c4.caption("Set conditions on *Describe experiment* to filter by replicates.")
        st.caption("Hide or show contaminants with the switch in the table's title row.")
    return query, q_max, int(min_peptides), int(min_replicates)
