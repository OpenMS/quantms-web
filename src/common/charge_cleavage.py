"""Precursor charge states and missed cleavages of the identified peptides, per file.

Both say whether the sample preparation and the acquisition went as planned.
A tryptic digest of a typical sample gives mostly doubly and triply charged
precursors; a run dominated by singly charged ones points to contamination
(polymers, detergents) or to a charge filter the instrument method did not
apply. Missed cleavages measure the digestion: most peptides should have none,
and a high share with one or more says the digestion was incomplete (too
little enzyme, too short, or a denaturation problem).

Reads the IDFilter idXMLs (target PSMs at the run's PSM FDR). The enzyme is
the one each file was searched with (its search parameters), so a run digested
with another protease counts its own cleavage sites. Streamlit-free, so tests
call it.
"""
from pathlib import Path

import pandas as pd
import pyopenms as poms

from src.common.qc_plots import _is_decoy
from src.common.results_helpers import extract_filename_from_idxml, load_idxml

DEFAULT_ENZYME = "Trypsin/P"
# Enzymes without cleavage sites: a missed cleavage means nothing there
NO_CLEAVAGE_SITES = {"unspecific cleavage", "no cleavage", "unspecific", "none", ""}

# The highest charge and missed cleavage counts that get a bar of their own;
# larger values are counted with them ("4+", "2+")
MAX_CHARGE = 4
MAX_MISSED = 2

PSM_COLUMNS = ["file", "sequence", "charge", "missed_cleavages"]
ALL_FILES = "All files"


def searched_enzyme(protein_ids) -> str:
    """The enzyme the file was searched with, from its search parameters."""
    for prot in protein_ids:
        try:
            name = prot.getSearchParameters().digestion_enzyme.getName()
        except (AttributeError, RuntimeError):
            continue
        # A file written without search parameters names "unknown_enzyme"
        if name and name != "unknown_enzyme":
            return name
    return DEFAULT_ENZYME


def _digestion(enzyme: str):
    """A ProteaseDigestion for ``enzyme``, None for an enzyme without cleavage sites."""
    if enzyme.lower() in NO_CLEAVAGE_SITES:
        return None
    digestion = poms.ProteaseDigestion()
    try:
        digestion.setEnzyme(enzyme)
    except RuntimeError:
        digestion.setEnzyme(DEFAULT_ENZYME)
    return digestion


def read_psms(idxml_path: Path, mzml_names=None) -> tuple[pd.DataFrame, str]:
    """The target PSMs (top hit of each spectrum) of one idXML: file, sequence,
    charge and missed cleavages, and the enzyme the file was searched with.
    Missed cleavages are empty for a search without an enzyme."""
    proteins, peptides = load_idxml(idxml_path)
    enzyme = searched_enzyme(proteins)
    digestion = _digestion(enzyme)
    file = extract_filename_from_idxml(Path(idxml_path), mzml_names)
    rows = []
    for pep in peptides:
        hits = pep.getHits()
        if not hits or _is_decoy(hits[0]):
            continue
        hit = hits[0]
        unmodified = hit.getSequence().toUnmodifiedString()
        rows.append({
            "file": file,
            "sequence": unmodified,
            "charge": hit.getCharge(),
            "missed_cleavages": digestion.countInternalCleavageSites(unmodified) if digestion else None,
        })
    return pd.DataFrame(rows, columns=PSM_COLUMNS), enzyme


def _charge_bin(charge: int) -> str:
    return f"{MAX_CHARGE}+" if charge >= MAX_CHARGE else str(int(charge))


def _missed_bin(missed: int) -> str:
    return f"{MAX_MISSED}+" if missed >= MAX_MISSED else str(int(missed))


CHARGE_BINS = [str(z) for z in range(1, MAX_CHARGE)] + [f"{MAX_CHARGE}+"]
MISSED_BINS = [str(m) for m in range(MAX_MISSED)] + [f"{MAX_MISSED}+"]


def _shares(values: pd.Series, bins: list[str]) -> pd.Series:
    """Share (%) of each bin, 0 for an empty bin."""
    counts = values.value_counts().reindex(bins, fill_value=0)
    total = counts.sum()
    return (counts / total * 100.0) if total else counts.astype(float)


def summary(psms: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Per-file table and the two distributions behind it.

    Returns ``(per_file, charges, missed)``. ``per_file`` has one row per file,
    with ``All files`` first (``file_id`` 0): target PSMs, the share of each
    charge state, the share of peptides by missed cleavages and the mean number
    of missed cleavages. ``charges`` and ``missed`` hold, per ``file_id`` and
    bin, the share (%) and the PSM count, for one bar plot each.
    Missed cleavages are counted over the distinct peptide sequences of a file,
    so a peptide identified in many spectra counts once.
    """
    psms = psms.copy()
    psms["charge_bin"] = psms["charge"].map(_charge_bin)
    has_missed = psms["missed_cleavages"].notna().any()
    files = sorted(psms["file"].unique())
    groups = [(0, ALL_FILES, psms)] + [(i + 1, f, psms[psms["file"] == f]) for i, f in enumerate(files)]
    rows, charge_rows, missed_rows = [], [], []
    for file_id, name, part in groups:
        charge = _shares(part["charge_bin"], CHARGE_BINS)
        counts = part["charge_bin"].value_counts()
        row = {"file_id": file_id, "File": name, "Target PSMs": len(part)}
        for bin_, share in charge.items():
            row[f"Charge {bin_} (%)"] = round(float(share), 1)
            charge_rows.append({
                "file_id": file_id, "bin": int(bin_.rstrip("+")), "label": f"{bin_}",
                "share": float(share), "psms": int(counts.get(bin_, 0)),
            })
        if has_missed:
            peptides = part.drop_duplicates(["file", "sequence"])
            missed_values = peptides["missed_cleavages"].dropna().astype(int)
            missed = _shares(missed_values.map(_missed_bin), MISSED_BINS)
            missed_counts = missed_values.map(_missed_bin).value_counts()
            for bin_, share in missed.items():
                row[f"{bin_} missed (%)"] = round(float(share), 1)
                missed_rows.append({
                    "file_id": file_id, "bin": int(bin_.rstrip("+")), "label": f"{bin_}",
                    "share": float(share), "peptides": int(missed_counts.get(bin_, 0)),
                })
            row["Mean missed"] = round(float(missed_values.mean()), 2) if len(missed_values) else None
        rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(charge_rows), pd.DataFrame(missed_rows)


# Sequential blues: low charge / no missed cleavage light, high dark
_SHADES = ["#c6dbef", "#6baed6", "#2a78d6", "#08306b"]


def stacked_shares(data: pd.DataFrame, per_file: pd.DataFrame, title: str, legend_title: str, count: str):
    """One stacked bar per file (All files first): the share of each bin, adding up to 100 %."""
    import plotly.graph_objects as go

    names = per_file.set_index("file_id")["File"]
    order = list(per_file["file_id"])
    fig = go.Figure()
    bins = data.drop_duplicates("bin").sort_values("bin")
    for color, (_, b) in zip(_SHADES[-len(bins):] if len(bins) < len(_SHADES) else _SHADES, bins.iterrows()):
        part = data[data["bin"] == b["bin"]].set_index("file_id").reindex(order)
        fig.add_trace(go.Bar(
            y=[names[i] for i in order], x=part["share"], name=b["label"], orientation="h",
            marker_color=color,
            customdata=part[count],
            hovertemplate=f"%{{y}}<br>{legend_title} {b['label']}: %{{x:.1f}} % (%{{customdata:,}} {count})<extra></extra>",
            text=[f"{v:.0f} %" if v >= 6 else "" for v in part["share"]], textposition="inside",
        ))
    fig.update_layout(
        barmode="stack", title=title, xaxis_title="Share (%)", xaxis_range=[0, 100],
        yaxis=dict(autorange="reversed"), legend=dict(title=legend_title, orientation="h", yanchor="bottom", y=1.02, x=0, traceorder="normal"),
        height=max(220, 34 * len(order) + 150), margin=dict(l=10),
    )
    return fig
