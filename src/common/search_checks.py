"""Checks of the search settings against the selected runs and database.

Streamlit-free. Before a search starts, Describe experiment compares what the
files say (how their MS2 spectra were fragmented and read out, how large the
protein database is) with the settings, and lists each mismatch as a notice
the user acts on or accepts. Nothing here stops a run.
"""
from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from src.common import instrument
from src.common.contaminants import CONTAMINANT_PREFIX

# Fragmentation, by the PSI-MS name of the activation of an MS2 spectrum
HCD, CID, ETD = "HCD", "CID", "ETD"
_ACTIVATIONS = {
    "beam-type collision-induced dissociation": HCD,
    "higher energy beam-type collision-induced dissociation": HCD,
    "collision-induced dissociation": CID,
    "trap-type collision-induced dissociation": CID,
    "electron transfer dissociation": ETD,
    "electron capture dissociation": ETD,
    "electron-transfer/higher-energy collision dissociation": ETD,
    "electron-transfer/collision-induced dissociation": ETD,
}
ACTIVATION_LABELS = {HCD: "HCD", CID: "CID", ETD: "ETD or EThcD"}
_ACTIVATION = re.compile(r"<activation>(.*?)</activation>", re.S)
_CV_NAME = re.compile(r'name="([^"]+)"')
_FILTER = re.compile(r'name="filter string" value="([^"]*)"')
_FILTER_ACTIVATION = re.compile(r"@(hcd|cid|etd|ethcd|etcid)", re.I)

# Sizes of a protein database (target entries)
LARGE_DATABASE = 150_000   # more than any single reference proteome
LARGE_OPEN_DATABASE = 40_000  # an open search compares every spectrum with far more peptides

# Fragment tolerances (Da) that fit an analyzer: ion trap MS2 needs about 0.5 Da
ION_TRAP_MIN_DA = 0.2
HIGH_RES_MAX_DA = 0.1
DECOY_PREFIXES = ("DECOY_", "decoy_", "rev_", "REV_", "XXX_", "reversed_")
_ORGANISM = re.compile(r"\bOS=(.+?)(?:\s+[A-Z]{2}=|$)")


@dataclass
class RunDetails:
    """What one run's file says about its MS2 spectra."""
    name: str
    model: str | None = None
    analyzer: str | None = None   # instrument preset of the MS2 readout
    activation: str | None = None  # HCD, CID or ETD, the most common in the file's head


@dataclass
class DatabaseSummary:
    """The proteins of a FASTA file."""
    name: str
    proteins: int = 0       # target entries, contaminants included
    decoys: int = 0
    contaminants: int = 0
    residues: int = 0
    organisms: dict = field(default_factory=dict)  # the five most common OS= names
    duplicates: int = 0
    readable: bool = True

    def text(self) -> str:
        if not self.readable:
            return f"**{self.name}** could not be read."
        if not self.proteins:
            return f"**{self.name}** holds no protein sequences."
        parts = [f"**{self.name}**: {self.proteins:,} proteins"]
        if self.organisms:
            names = [f"*{o}* ({n:,})" for o, n in self.organisms.items()]
            top = names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]
            parts.append(f"from {top}")
        extras = []
        if self.contaminants:
            extras.append(f"{self.contaminants:,} contaminants")
        if self.decoys:
            extras.append(f"{self.decoys:,} decoys already in the file")
        text = " ".join(parts)
        if extras:
            text += f" ({', '.join(extras)})"
        return text + "."


@dataclass
class Notice:
    """One mismatch: what is wrong, and the preset that fixes it, if one does."""
    id: str
    level: str  # "warning" or "info"
    text: str
    preset: str | None = None
    fingerprint: str = ""


def _activation_of(block: str) -> str | None:
    for name in _CV_NAME.findall(block):
        if name in _ACTIVATIONS:
            return _ACTIVATIONS[name]
    return None


def read_activation(head: str) -> str | None:
    """The fragmentation most MS2 spectra in an mzML head were made with."""
    counts: Counter = Counter()
    for block in _ACTIVATION.findall(head):
        kind = _activation_of(block)
        if kind:
            counts[kind] += 1
    if not counts:
        for filt in _FILTER.findall(head):
            match = _FILTER_ACTIVATION.search(filt)
            if match:
                word = match.group(1).lower()
                counts[ETD if word.startswith("et") else HCD if word == "hcd" else CID] += 1
    return counts.most_common(1)[0][0] if counts else None


def read_run(path: str | Path, converted_dir: str | Path | None = None) -> RunDetails | None:
    """Details of one selected run, or None when its file cannot be read yet."""
    path = Path(path)
    if path.suffix.lower() == ".d":
        return RunDetails(path.name, "Bruker timsTOF", instrument.TIMSTOF, None)
    readable = path
    if path.suffix.lower() == ".raw":
        if converted_dir is None:
            return None
        readable = Path(converted_dir, f"{path.stem}.mzML")
    if not readable.is_file():
        return None
    try:
        with open(readable, "rb") as f:
            head = f.read(instrument.HEAD_BYTES).decode("utf-8", errors="ignore")
    except OSError:
        return None
    run = instrument.read_mzml(readable)
    return RunDetails(path.name, run.model, run.preset, read_activation(head))


def summarize_database(path: str | Path) -> DatabaseSummary:
    """Count the target, decoy and contaminant entries of a FASTA file."""
    path = Path(path)
    summary = DatabaseSummary(path.name)
    organisms: Counter = Counter()
    seen: set = set()
    try:
        with open(path, "r", errors="ignore") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip()
                    accession = header.split(None, 1)[0] if header else ""
                    if accession in seen:
                        summary.duplicates += 1
                    seen.add(accession)
                    if accession.startswith(DECOY_PREFIXES):
                        summary.decoys += 1
                        continue
                    summary.proteins += 1
                    if accession.startswith(CONTAMINANT_PREFIX):
                        summary.contaminants += 1
                    else:
                        match = _ORGANISM.search(header)
                        if match:
                            organisms[match.group(1).split(" (", 1)[0].strip()] += 1
                elif not line.startswith(";"):
                    summary.residues += len(line.strip())
    except OSError:
        summary.readable = False
        return summary
    summary.organisms = dict(organisms.most_common(5))
    return summary


def run_summary(runs: list) -> str:
    """The instruments and fragmentation of the readable runs, in one line."""
    groups: Counter = Counter()
    for run in runs:
        what = []
        if run.activation:
            what.append(ACTIVATION_LABELS[run.activation])
        if run.analyzer == instrument.ION_TRAP:
            what.append("MS2 in the ion trap")
        elif run.analyzer == instrument.ORBITRAP:
            what.append("MS2 in the Orbitrap")
        name = run.model or "unknown instrument"
        groups[f"{name} ({', '.join(what)})" if what else name] += 1
    return "; ".join(f"{n} × {w}" for w, n in groups.most_common())


def _fragment_da(values: dict) -> tuple[float, str]:
    tolerance = float(values.get("Search:fragment:mass_tolerance", 20.0) or 0.0)
    unit = values.get("Search:fragment:mass_tolerance_unit", "ppm")
    return tolerance, unit


def _plural(n: int, word: str) -> str:
    return f"{n} {word}{'' if n == 1 else 's'}"


def _share(part: int, whole: int) -> str:
    """'The run', 'Both runs', 'All 6 runs' or '2 of 6 runs', to start a sentence."""
    if part == whole:
        return "The run" if whole == 1 else "Both runs" if whole == 2 else f"All {whole} runs"
    return f"{part} of {whole} runs"


def check(runs: list, values: dict, database: DatabaseSummary | None, open_search: bool = False) -> list:
    """The notices for the readable runs, the search settings in effect and the database."""
    notices = []
    tolerance, unit = _fragment_da(values)
    shown = f"{tolerance:g} {unit}"

    trap = [r for r in runs if r.analyzer == instrument.ION_TRAP]
    if trap and (unit == "ppm" or tolerance < ION_TRAP_MIN_DA):
        notices.append(Notice(
            "fragment-ion-trap", "warning",
            f"{_share(len(trap), len(runs))} read their MS2 spectra in the ion trap, whose "
            f"fragment masses are only accurate to about 0.5 Da, but the fragment tolerance is {shown}. "
            "Most fragments will fall outside it and few peptides will be identified.",
            instrument.ION_TRAP, f"{len(trap)}|{shown}",
        ))
    high_res = [r for r in runs if r.analyzer in (instrument.ORBITRAP, instrument.TIMSTOF)]
    if high_res and not trap and unit == "Da" and tolerance > HIGH_RES_MAX_DA:
        preset = high_res[0].analyzer
        notices.append(Notice(
            "fragment-high-res", "info",
            f"{_share(len(high_res), len(runs))} measured their fragments at high resolution, "
            f"but the fragment tolerance is {shown}, wide enough for a peak of another fragment to match. "
            "A tolerance in ppm finds more peptides at the same error rate.",
            preset, f"{len(high_res)}|{shown}",
        ))

    etd = [r for r in runs if r.activation == ETD]
    if etd:
        notices.append(Notice(
            "fragmentation-etd", "warning",
            f"{_share(len(etd), len(runs))} fragmented their peptides with ETD or EThcD, "
            "which breaks them into c and z ions. The search looks for the b and y ions of HCD and CID, "
            "so most of these spectra will stay unidentified.",
            None, ",".join(sorted(r.name for r in etd)),
        ))
    kinds = {r.activation for r in runs if r.activation}
    if len(kinds) > 1 and not etd:
        notices.append(Notice(
            "fragmentation-mixed", "info",
            "The selected runs were fragmented in different ways ("
            + ", ".join(sorted(ACTIVATION_LABELS[k] for k in kinds))
            + "), and are searched with the same settings. Check that they belong to one experiment.",
            None, ",".join(sorted(kinds)),
        ))

    if database is not None:
        if not database.readable or not database.proteins:
            notices.append(Notice(
                "database-empty", "warning",
                f"{database.text()} Choose a FASTA file of your organism's proteins with **📂 Browse**.",
                None, database.name,
            ))
        else:
            limit = LARGE_OPEN_DATABASE if open_search else LARGE_DATABASE
            if database.proteins > limit:
                what = "An open search" if open_search else "The search"
                notices.append(Notice(
                    "database-large", "warning",
                    f"**{database.name}** holds {database.proteins:,} proteins, more than a reference "
                    "proteome (4,000 to 25,000 for most organisms, about 20,000 for human). "
                    f"{what} against it takes much longer, and in a larger search space fewer matches pass "
                    "the FDR filter. A reference proteome of your organism is under "
                    "**📂 Browse** → *Import from UniProt proteome*.",
                    None, f"{database.name}|{database.proteins}|{open_search}",
                ))
            if database.duplicates:
                notices.append(Notice(
                    "database-duplicates", "info",
                    f"**{database.name}** lists {_plural(database.duplicates, 'accession')} more than once. "
                    "Often two databases were joined that share entries; check that it is the file you meant.",
                    None, f"{database.name}|{database.duplicates}",
                ))
    return notices
