"""Precursor mass shifts of identified PSMs, and the Unimod modifications they match.

In an open search (precursor tolerance of hundreds of Da) a PSM's precursor may
carry a modification the search did not consider. The shift between observed
and theoretical precursor mass then shows which ones: shifts recurring across
many PSMs form peaks, and each peak is matched to the Unimod modifications of
that mass (pyOpenMS ModificationsDB). In a closed search the same view shows
the 0 Da peak and the isotope picks.

Reads only the IDFilter idXMLs. It makes no Streamlit calls, so the page caches
it and tests call it.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from pyopenms import ModificationsDB, ResidueModification

from src.common.results_helpers import extract_filename_from_idxml, load_idxml

PROTON = 1.007276466812
C13_DELTA = 1.0033548378

# Width of a histogram bin, and the gap that always separates two shifts;
# closer peaks (deamidation +0.984 Da, a 13C pick +1.003 Da) are split at their valley
BIN_WIDTH = 0.005
CLUSTER_GAP = 0.02
# A shift needs this many PSMs to be reported
MIN_PSMS = 3
# Unimod matches must be this close to the observed shift
MATCH_TOLERANCE = 0.02

# Shown first when several modifications share a mass
COMMON_MODIFICATIONS = [
    "Oxidation", "Phospho", "Carbamidomethyl", "Acetyl", "Deamidated", "Methyl",
    "Dimethyl", "Trimethyl", "GG", "Amidated", "Ammonia-loss", "Gln->pyro-Glu",
    "Glu->pyro-Glu", "Carbamyl", "Formyl", "Sodium", "Potassium", "Dehydrated",
    "Dioxidation", "Trioxidation", "Nitro", "TMT6plex", "TMTpro", "iTRAQ4plex",
]

# Shifts that are not modifications: a monoisotopic peak missed by n isotopes
ISOTOPE_LABELS = {1: "13C isotope peak", 2: "2 × 13C isotope peak", 3: "3 × 13C isotope peak"}

SHIFT_COLUMNS = [
    "shift_id", "Observed Δ (Da)", "PSMs", "Interpretation", "Unimod Δ (Da)",
    "Deviation (mDa)", "Residues", "Other matches",
]


def _is_decoy(hit) -> bool:
    if not hit.metaValueExists("target_decoy"):
        return False
    value = hit.getMetaValue("target_decoy")
    if isinstance(value, bytes):
        value = value.decode()
    return str(value).startswith("decoy")


def psm_mass_shifts(idxml_paths, mzml_names=None) -> pd.DataFrame:
    """Target top hits with their precursor mass shift (observed − theoretical, in Da)."""
    rows = []
    for path in idxml_paths:
        source = extract_filename_from_idxml(Path(path), mzml_names)
        _, peptides = load_idxml(path)
        for pep in peptides:
            hits = pep.getHits()
            if not hits or _is_decoy(hits[0]):
                continue
            hit = hits[0]
            charge = hit.getCharge()
            if charge <= 0:
                continue
            sequence = hit.getSequence()
            observed = pep.getMZ() * charge - charge * PROTON
            rows.append({
                "file": source,
                "sequence": sequence.toString(),
                "charge": charge,
                "mz": pep.getMZ(),
                "rt": pep.getRT(),
                "score": hit.getScore(),
                "delta": observed - sequence.getMonoWeight(),
            })
    return pd.DataFrame(rows, columns=["file", "sequence", "charge", "mz", "rt", "score", "delta"])


def _split_points(values: np.ndarray) -> list[float]:
    """Mass shifts at which one gap-free group of sorted shifts holds two peaks.

    The group is histogrammed in 1 mDa bins and smoothed; between two maxima it is
    split at the lowest point when that valley is below half the smaller maximum,
    so two overlapping peaks (deamidation and a 13C pick) separate while one
    noisy peak does not.
    """
    step = 0.001
    first = values.min()
    counts = np.bincount(((values - first) / step).astype(int))
    kernel = np.exp(-0.5 * (np.arange(-6, 7) / 2.0) ** 2)
    smooth = np.convolve(counts, kernel / kernel.sum(), mode="same")
    maxima = [
        i for i in range(len(smooth))
        if smooth[i] > 0
        and (i == 0 or smooth[i] > smooth[i - 1])
        and (i == len(smooth) - 1 or smooth[i] >= smooth[i + 1])
    ]
    splits = []
    kept = maxima[:1]
    for i in maxima[1:]:
        left = kept[-1]
        valley = left + int(np.argmin(smooth[left:i + 1]))
        if smooth[valley] < 0.5 * min(smooth[left], smooth[i]):
            splits.append(first + (valley + 0.5) * step)
            kept.append(i)
        elif smooth[i] > smooth[left]:
            kept[-1] = i
    return splits


def _clusters(deltas: np.ndarray) -> list[np.ndarray]:
    """Indices of ``deltas`` grouped into shifts.

    Sorted shifts more than CLUSTER_GAP apart start a new group, and a group with
    two peaks is split between them.
    """
    if len(deltas) == 0:
        return []
    order = np.argsort(deltas)
    breaks = np.where(np.diff(deltas[order]) > CLUSTER_GAP)[0] + 1
    clusters = []
    for group in np.split(order, breaks):
        values = deltas[group]
        cuts = _split_points(values) if len(group) >= 2 * MIN_PSMS else []
        position = np.searchsorted(cuts, values)
        clusters.extend(group[position == k] for k in range(len(cuts) + 1))
    return [c for c in clusters if len(c)]


def unimod_matches(shift: float, tolerance: float = MATCH_TOLERANCE) -> list[dict]:
    """Unimod modifications within ``tolerance`` of ``shift``, one entry per name.

    Each entry has the name, its Δ mass, the deviation of the observed shift from
    it and the residues it is defined on; common modifications come first, then
    the closest.
    """
    full_ids = []
    ModificationsDB().searchModificationsByDiffMonoMass(
        full_ids, float(shift), float(tolerance), "", ResidueModification.TermSpecificity.ANYWHERE
    )
    by_name = {}
    db = ModificationsDB()
    for full_id in full_ids:
        full_id = full_id.decode() if isinstance(full_id, bytes) else str(full_id)
        mod = db.getModification(full_id)
        name = mod.getId()
        # Run-time placeholders and PSI-MOD duplicates of Unimod entries
        if not name or name.startswith(("CUSTOM", "[", "MOD:")):
            continue
        entry = by_name.setdefault(name, {
            "name": name,
            "delta": mod.getDiffMonoMass(),
            "residues": set(),
        })
        origin = mod.getOrigin()
        term = mod.getTermSpecificity()
        if term == ResidueModification.TermSpecificity.N_TERM:
            entry["residues"].add("N-term")
        elif term == ResidueModification.TermSpecificity.C_TERM:
            entry["residues"].add("C-term")
        elif origin and origin != "X":
            entry["residues"].add(origin)
    matches = []
    for entry in by_name.values():
        entry["deviation"] = shift - entry["delta"]
        entry["residues"] = "".join(sorted(r for r in entry["residues"] if len(r) == 1)) + "".join(
            f" {r}" for r in sorted(r for r in entry["residues"] if len(r) > 1)
        )
        matches.append(entry)
    rank = {name: i for i, name in enumerate(COMMON_MODIFICATIONS)}
    matches.sort(key=lambda m: (rank.get(m["name"], len(rank)), abs(m["deviation"]), m["name"]))
    return matches


def _interpret(shift: float) -> dict:
    """Interpretation of one observed shift: unmodified, isotope pick, or Unimod matches.

    A shift close to both an isotope pick and a modification (deamidation, +0.984 Da,
    is 19 mDa from one 13C) goes to whichever is closer; the other is listed too.
    """
    if abs(shift) <= MATCH_TOLERANCE:
        return {"Interpretation": "Unmodified", "Unimod Δ (Da)": 0.0,
                "Deviation (mDa)": round(shift * 1000, 1), "Residues": "", "Other matches": ""}
    candidates = [
        {"name": label, "delta": n * C13_DELTA, "deviation": shift - n * C13_DELTA, "residues": ""}
        for n, label in ISOTOPE_LABELS.items()
        if abs(shift - n * C13_DELTA) <= MATCH_TOLERANCE
    ]
    matches = unimod_matches(shift)
    if candidates and (not matches or abs(candidates[0]["deviation"]) < abs(matches[0]["deviation"])):
        matches = candidates + matches
    else:
        matches = matches + candidates
    if not matches:
        return {"Interpretation": "", "Unimod Δ (Da)": None, "Deviation (mDa)": None,
                "Residues": "", "Other matches": ""}
    best = matches[0]
    return {
        "Interpretation": best["name"],
        "Unimod Δ (Da)": round(best["delta"], 4),
        "Deviation (mDa)": round(best["deviation"] * 1000, 1),
        "Residues": best["residues"].strip(),
        "Other matches": ", ".join(m["name"] for m in matches[1:6]),
    }


def mass_shift_summary(psms: pd.DataFrame, min_psms: int = MIN_PSMS) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Detected shifts, the histogram and the PSMs, linked by ``shift_id``.

    Returns:
        shifts: SHIFT_COLUMNS, one row per shift with at least ``min_psms`` PSMs,
            most PSMs first.
        histogram: one row per non-empty BIN_WIDTH bin (``delta``, ``psms``),
            with the tallest bin of each shift marked (``peak``) and labelled.
        psms: the input PSMs with the ``shift_id`` they belong to (-1 for none).
    """
    psms = psms.copy()
    psms["shift_id"] = -1
    deltas = psms["delta"].to_numpy(dtype=float)

    shift_rows = []
    for members in _clusters(deltas):
        if len(members) < min_psms:
            continue
        observed = float(np.median(deltas[members]))
        shift_rows.append({"members": members, "Observed Δ (Da)": round(observed, 4), "PSMs": len(members),
                           **_interpret(observed)})
    shift_rows.sort(key=lambda r: (-r["PSMs"], r["Observed Δ (Da)"]))
    for shift_id, row in enumerate(shift_rows):
        row["shift_id"] = shift_id
        psms.iloc[row.pop("members"), psms.columns.get_loc("shift_id")] = shift_id
    shifts = pd.DataFrame(shift_rows, columns=SHIFT_COLUMNS)
    for col in ("Observed Δ (Da)", "Unimod Δ (Da)", "Deviation (mDa)"):
        shifts[col] = pd.to_numeric(shifts[col], errors="coerce")

    psms["bin"] = np.round(deltas / BIN_WIDTH).astype(int) if len(psms) else pd.Series(dtype=int)
    histogram = (
        psms.groupby("bin")
        .agg(psms=("delta", "size"), shift_id=("shift_id", "max"))
        .reset_index()
    )
    histogram["delta"] = histogram["bin"] * BIN_WIDTH
    histogram["peak"] = False
    histogram["label"] = ""
    labels = dict(zip(shifts["shift_id"], shifts["Interpretation"]))
    observed = dict(zip(shifts["shift_id"], shifts["Observed Δ (Da)"]))
    for shift_id, group in histogram[histogram["shift_id"] >= 0].groupby("shift_id"):
        top = group["psms"].idxmax()
        histogram.loc[top, "peak"] = True
        label = labels.get(shift_id) or "?"
        shift = observed[shift_id] if abs(observed[shift_id]) >= 0.0005 else 0.0
        histogram.loc[top, "label"] = f"{shift:+.3f} {label}"
    return shifts, histogram[["bin", "delta", "psms", "shift_id", "peak", "label"]], psms.drop(columns="bin")
