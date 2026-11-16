"""Fragment-level evidence for one PSM: its matched fragments and the b/y ion table
(the mass error plot is the library's ``identification/fragment_errors``).

Both read the run's ``<prefix>_fragment_ions.tsv`` export (written by the
workflow with the same ion types and tolerance as the SequenceView), filtered
to the selected PSM.
"""
import html
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from pyopenms import AASequence, Residue

from src.viz.identification.fragment_errors import ION_COLORS

# The ion table marks matches in the colours the fragment error plot gives them
B_COLOR = ION_COLORS["b"]
Y_COLOR = ION_COLORS["y"]


def load_psm_fragments(tsv_path: Path, sequence_id: int) -> pd.DataFrame:
    """The fragment ions matched in one PSM's spectrum, read lazily from the run's export."""
    return (
        pl.scan_csv(tsv_path, separator="\t")
        .filter(pl.col("sequence_id") == sequence_id)
        .collect()
        .to_pandas()
    )


def _is_plain(ion: str, ion_type: str, number: int) -> bool:
    """Whether an ion label is the plain fragment, not a neutral-loss one (``b4`` vs ``b4-H2O``)."""
    return ion == f"{ion_type}{number}"


def ion_table(sequence: str, precursor_charge: int, fragments: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Theoretical b and y m/z per position, and which of them a peak matched.

    One row per residue; b ions count from the N-terminus and y ions from the
    C-terminus, at charges 1 up to ``precursor_charge - 1`` (at most 2). Only
    plain ions count as matches, not their neutral losses.

    Returns:
        (values, matched): the m/z table and a boolean table of the same shape.
    """
    peptide = AASequence.fromString(sequence)
    n = peptide.size()
    charges = list(range(1, max(1, min(int(precursor_charge) - 1, 2)) + 1))
    plain = fragments[[
        _is_plain(i, t, k) for i, t, k in zip(fragments["ion"], fragments["ion_type"], fragments["ion_number"])
    ]] if len(fragments) else fragments
    matched_keys = set(zip(plain["ion_type"], plain["ion_number"], plain["charge"])) if len(plain) else set()

    def column(ion_type, z):
        return f"{ion_type}{'+' * z}"

    values = {"#": list(range(1, n + 1)), "residue": [peptide.getResidue(i).getOneLetterCode() for i in range(n)]}
    matched = {"#": [False] * n, "residue": [False] * n}
    for z in charges:
        b_vals, b_hit, y_vals, y_hit = [], [], [], []
        for i in range(1, n + 1):
            # b_i covers residues 1..i (no b_n), y_i the last i residues; row i shows b_i and y_(n-i+1)
            if i < n:
                b_vals.append(peptide.getPrefix(i).getMonoWeight(Residue.ResidueType.BIon, z) / z)
                b_hit.append(("b", i, z) in matched_keys)
            else:
                b_vals.append(np.nan)
                b_hit.append(False)
            k = n - i + 1
            if k < n:
                y_vals.append(peptide.getSuffix(k).getMonoWeight(Residue.ResidueType.YIon, z) / z)
                y_hit.append(("y", k, z) in matched_keys)
            else:
                y_vals.append(np.nan)
                y_hit.append(False)
        values[column("b", z)], matched[column("b", z)] = b_vals, b_hit
        values[column("y", z)], matched[column("y", z)] = y_vals, y_hit
    values["#y"] = [n - i + 1 for i in range(1, n + 1)]
    matched["#y"] = [False] * n
    order = ["#"] + [column("b", z) for z in charges] + ["residue"] + [column("y", z) for z in charges] + ["#y"]
    return pd.DataFrame(values)[order], pd.DataFrame(matched)[order]


def ion_table_html(values: pd.DataFrame, matched: pd.DataFrame) -> str:
    """The ion table as HTML, matched b cells in blue and y cells in orange."""
    head = "".join(f"<th>{html.escape('#' if c == '#y' else c)}</th>" for c in values.columns)
    rows = []
    for r in range(len(values)):
        cells = []
        for c in values.columns:
            v = values.iat[r, values.columns.get_loc(c)]
            text = "" if isinstance(v, float) and np.isnan(v) else (f"{v:.4f}" if isinstance(v, float) else html.escape(str(v)))
            cls = ""
            if matched.iat[r, matched.columns.get_loc(c)]:
                cls = ' class="hit-b"' if c.startswith("b") else ' class="hit-y"'
            elif c == "residue":
                cls = ' class="res"'
            cells.append(f"<td{cls}>{text}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return (
        "<style>.iontab{border-collapse:collapse;font-family:monospace;font-size:13px}"
        ".iontab th,.iontab td{padding:2px 10px;text-align:right;border-bottom:1px solid #e4e4e7}"
        ".iontab .res{text-align:center;font-weight:700}"
        f".iontab .hit-b{{background:{B_COLOR};color:white}}.iontab .hit-y{{background:{Y_COLOR};color:white}}</style>"
        ".iontab caption{caption-side:top;text-align:left;font-family:inherit;font-weight:600;padding:4px 0}"
        f'<table class="iontab"><caption>b and y ions (theoretical m/z, matched ions highlighted)</caption><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table>'
    )
