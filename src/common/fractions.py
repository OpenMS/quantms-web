"""Samples of a fractionated run: which input files are fractions of one sample.

ProteomicsLFQ gets the experimental design when a sample was split into
fractions, but its MSstats table still has one ``Reference`` (input file) per
row, so read per file the results show every fraction as its own sample. The
results pages use :func:`fraction_samples` to put the fractions of a sample
back together.

A run is fractionated when the workflow passed a design to ProteomicsLFQ, which
it saves as ``quant_results/experimental-design.tsv``. Files that share a
condition and biological replicate are one sample; a second set of fractions of
the same sample (a technical replicate) is a sample of its own.
"""
from pathlib import Path
from typing import Dict, Mapping, Optional

import pandas as pd

from src.workflow.ExperimentalDesign import DESIGN_FILE, read_design


def design_samples(design: pd.DataFrame, params: Mapping) -> Dict[str, dict]:
    """Each file's sample name and condition, for a design with fractions.

    The condition comes from the run's ``mzML-group-<file>`` parameter of the
    sample's first fraction (that is where the results pages read conditions
    from), else from the design. Samples are named ``<condition>_<replicate>``,
    with ``_t<n>`` added for technical replicates.
    """
    design = design.reset_index(drop=True).copy()
    design["technical"] = design.groupby(["condition", "replicate", "fraction"]).cumcount() + 1
    out = {}
    for (cond, rep, tech), grp in design.groupby(["condition", "replicate", "technical"], sort=False):
        first = grp.sort_values("fraction").iloc[0]["file"]
        condition = str(params.get(f"mzML-group-{first}") or cond or "").strip()
        name = f"{condition or 'sample'}_{rep}"
        if design[(design["condition"] == cond) & (design["replicate"] == rep)]["technical"].max() > 1:
            name += f"_t{tech}"
        for f in grp["file"]:
            out[f] = {"sample": name, "condition": condition}
    return out


def fraction_samples(results_dir, params: Mapping) -> Optional[Dict[str, dict]]:
    """``{file name: {"sample", "condition"}}`` for a fractionated run, else None.

    ``results_dir`` is the run's results folder (holding quant_results/).
    """
    path = Path(results_dir, "quant_results", DESIGN_FILE)
    if not path.is_file():
        return None
    try:
        design = read_design(path)
    except (OSError, ValueError):
        return None
    if design.empty or not (design["fraction"] > 1).any():
        return None
    return design_samples(design, params)


def sample_of(reference: str, samples: Optional[Mapping[str, dict]]) -> str:
    """The sample an MSstats ``Reference`` belongs to (the file's stem without fractions)."""
    name = Path(str(reference)).name
    if samples and name in samples:
        return samples[name]["sample"]
    return name[:-5] if name.endswith(".mzML") else Path(name).stem
