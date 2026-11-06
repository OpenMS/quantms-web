"""Which instrument preset fits the selected runs, read from their files.

Streamlit-free, so it can be tested on its own. The instrument is read from
the head of an mzML file: the model and mass analyzers in its
``instrumentConfigurationList`` and the scan filters of its first MS2
spectra (``FTMS`` = Orbitrap, ``ITMS`` = ion trap). A Bruker ``.d`` folder is
a timsTOF run; a Thermo ``.raw`` file is read from its converted mzML once
it has been converted, and is not counted before that.
"""
from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

ORBITRAP = "Orbitrap HCD"
TIMSTOF = "timsTOF"
ION_TRAP = "Ion Trap CID"
PRESETS = (ORBITRAP, TIMSTOF, ION_TRAP)

# How much of an mzML file is read: the header and the first spectra.
HEAD_BYTES = 4 * 1024 * 1024

_MODEL_ACCESSIONS = re.compile(r'accession="MS:(\d+)" name="([^"]+)"')
_FILTER = re.compile(r'name="filter string" value="([^"]*)"')
_CONFIG = re.compile(r"<instrumentConfigurationList.*?</instrumentConfigurationList>", re.S)
_ANALYZERS = {
    "orbitrap": ORBITRAP,
    "radial ejection linear ion trap": ION_TRAP,
    "axial ejection linear ion trap": ION_TRAP,
    "quadrupole ion trap": ION_TRAP,
    "time-of-flight": TIMSTOF,
}


@dataclass
class RunInstrument:
    """What one file says about the instrument that measured it."""
    model: str | None = None
    preset: str | None = None
    why: str = ""


@dataclass
class Recommendation:
    """The preset most selected runs fit, and why."""
    preset: str | None
    runs: int = 0
    unreadable: int = 0
    models: list = field(default_factory=list)
    counts: dict = field(default_factory=dict)

    def reason(self) -> str:
        if self.preset is None:
            if self.counts:
                parts = ", ".join(f"{n} {'run' if n == 1 else 'runs'} {p}" for p, n in self.counts.items())
                return f"The selected runs come from different instruments ({parts}): check the preset."
            return ""
        n = self.counts.get(self.preset, 0)
        what = f"{', '.join(self.models)}" if self.models else "the instrument"
        return (f"{n} of {self.runs} {'run was' if self.runs == 1 else 'runs were'} measured on "
                f"{what}, which the **{self.preset}** preset fits.")


def _model_name(config: str) -> str | None:
    """The instrument model: a term of the configuration outside its components."""
    outside = re.sub(r"<componentList.*?</componentList>", "", config, flags=re.S)
    for _, name in _MODEL_ACCESSIONS.findall(outside):
        if not any(s in name.lower() for s in ("serial", "customization", "instrument model", "software")):
            return name
    return None


def read_mzml(path: str | Path) -> RunInstrument:
    """Read the instrument of an mzML file from its head."""
    try:
        with open(path, "rb") as f:
            head = f.read(HEAD_BYTES).decode("utf-8", errors="ignore")
    except OSError:
        return RunInstrument()
    config_match = _CONFIG.search(head)
    config = config_match.group(0) if config_match else ""
    model = _model_name(config)

    # The analyzer of the MS2 scans decides, when the scan filters name it
    for filt in _FILTER.findall(head):
        if " ms2 " not in f" {filt} ":
            continue
        if filt.startswith("ITMS"):
            return RunInstrument(model, ION_TRAP, "MS2 in the ion trap")
        if filt.startswith("FTMS"):
            return RunInstrument(model, ORBITRAP, "MS2 in the Orbitrap")

    lowered = (model or "").lower()
    if "tims" in lowered:
        return RunInstrument(model, TIMSTOF, "timsTOF")
    analyzers = [p for name, p in _ANALYZERS.items() if f'name="{name}"' in config]
    if ORBITRAP in analyzers:
        return RunInstrument(model, ORBITRAP, "Orbitrap analyzer")
    if analyzers:
        return RunInstrument(model, analyzers[0], "mass analyzer")
    if "orbitrap" in lowered or "exactive" in lowered:
        return RunInstrument(model, ORBITRAP, "Orbitrap model")
    return RunInstrument(model)


def read_run(path: str | Path, converted_dir: str | Path | None = None) -> RunInstrument | None:
    """The instrument of one selected input, or None when it cannot be read yet."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".d":
        return RunInstrument("Bruker timsTOF", TIMSTOF, "Bruker .d folder")
    if suffix == ".raw":
        if converted_dir is None:
            return None
        converted = Path(converted_dir, f"{path.stem}.mzML")
        if not converted.is_file():
            return None
        path = converted
    run = read_mzml(path)
    return run if run.preset else None


def recommend(paths: list, converted_dir: str | Path | None = None) -> Recommendation:
    """The preset that fits all readable runs, or None when they disagree or none is readable."""
    counts: Counter = Counter()
    models: list = []
    unreadable = 0
    for p in paths:
        run = read_run(p, converted_dir)
        if run is None:
            unreadable += 1
            continue
        counts[run.preset] += 1
        if run.model and run.model not in models:
            models.append(run.model)
    runs = sum(counts.values())
    if not counts:
        return Recommendation(None, 0, unreadable)
    if len(counts) > 1:
        return Recommendation(None, runs, unreadable, models, dict(counts.most_common()))
    preset = next(iter(counts))
    return Recommendation(preset, runs, unreadable, models, dict(counts))
