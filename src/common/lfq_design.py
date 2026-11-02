"""Keep the LFQ experimental design and the per-file group parameters in step.

The results pages read each file's group from the ``mzML-group-<file>``
parameters, so those stay the source of the conditions. The design table adds
the biological replicate and fraction of each file and is saved as an OpenMS
design.tsv, which ProteomicsLFQ gets to combine the fractions of a sample.
"""
from typing import Callable, Dict, Mapping

import pandas as pd


def seed_design(
    design: pd.DataFrame, groups: Mapping[str, str], renumber: bool = False
) -> pd.DataFrame:
    """The design with each file's condition taken from ``groups``.

    With ``renumber`` (no design saved yet), the files of each condition are
    numbered as biological replicates 1, 2, 3, ... in their order.
    """
    out = design.copy()
    out["condition"] = [str(groups.get(f, "") or "") for f in out["file"]]
    if renumber:
        out["replicate"] = out.groupby("condition").cumcount() + 1
    return out


def groups_from_design(design: pd.DataFrame) -> Dict[str, str]:
    """Each file's condition, as stored in the ``mzML-group-<file>`` parameters."""
    return dict(zip(design["file"], design["condition"].fillna("").astype(str)))


def is_fractionated(design: pd.DataFrame) -> bool:
    """Whether any sample was split into fractions, which ProteomicsLFQ must combine.

    Without fractions every file is its own sample, which is what ProteomicsLFQ
    assumes when it gets no design, so the design is only passed when this is true.
    """
    return bool(not design.empty and (design["fraction"] > 1).any())


def initial_groups(files, params: Mapping, default: Callable[[str], str]) -> Dict[str, str]:
    """Each file's group from the parameters, or ``default(file)`` if it has none yet."""
    return {
        f: params[f"mzML-group-{f}"] if f"mzML-group-{f}" in params else default(f)
        for f in files
    }
