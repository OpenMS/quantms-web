"""Experiment-wide exports: every file's filtered PSMs in one idXML, and the design.tsv.

A run keeps one idXML per input file. These two files describe the whole
experiment at once, for tools that want one identification file and the design
that goes with it (OpenMS, MSstats, quantms):

- ``combined/all_psms_filtered.idXML``: IDMerger's merge of the per-file
  idXMLs after the PSM FDR filter. Each identification keeps its source file
  in the ``file_origin`` meta value.
- ``combined/design.tsv``: the OpenMS experimental design (two-table format)
  with bare file names, so it is valid on any machine.

This module has no Streamlit dependency, so ``execution()`` can use it in a
queue worker.
"""
import json
import re
import subprocess
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional

import pandas as pd

from src.common.isobaric import ISOBARIC_CHANNELS
from src.workflow.ExperimentalDesign import _file_name, align_design, read_design, to_openms_design

COMBINED_DIR = "combined"
MERGED_IDXML = "all_psms_filtered.idXML"
DESIGN_TSV = "design.tsv"


def combined_paths(results_dir) -> Dict[str, Path]:
    """Where a run's merged idXML and design.tsv go."""
    folder = Path(results_dir, COMBINED_DIR)
    return {"idxml": folder / MERGED_IDXML, "design": folder / DESIGN_TSV}


def lfq_design_tsv(design: pd.DataFrame, params: Mapping) -> str:
    """The label-free design as design.tsv, conditions taken from the run's parameters.

    The ``mzML-group-<file>`` parameters are where the results pages read the
    conditions from, so they win over the condition stored in ``design``.
    Spectra_Filepath holds bare file names.
    """
    design = design.copy()
    design["condition"] = [
        str(params.get(f"mzML-group-{f}", c) or "") for f, c in zip(design["file"], design["condition"])
    ]
    return to_openms_design(design)


def isobaric_design_tsv(files: List[str], params: Mapping, label: str) -> str:
    """An isobaric design.tsv: one row per file and channel, conditions from ``TMT-group-sample<n>``.

    Every file has the same channel layout, so channel n of every file is in
    the condition the user gave channel n. Each file and channel is its own
    sample (biological replicate). Channels marked "skip" or left empty keep
    the condition "NA".
    """
    m = re.search(r"\d+", label)
    plex = int(m.group()) if m else len(ISOBARIC_CHANNELS.get(label, []))
    names = [_file_name(str(f)) for f in files]
    file_rows, sample_rows = [], []
    for fi, name in enumerate(names, start=1):
        for ch in range(1, plex + 1):
            sample = (fi - 1) * plex + ch
            group = str(params.get(f"TMT-group-sample{ch}", "") or "").strip()
            if not group or group.lower() == "skip":
                group = "NA"
            file_rows.append({"Fraction_Group": fi, "Fraction": 1, "Spectra_Filepath": name, "Label": ch, "Sample": sample})
            sample_rows.append({"Sample": sample, "MSstats_Condition": group, "MSstats_BioReplicate": sample})
    return (
        pd.DataFrame(file_rows, columns=["Fraction_Group", "Fraction", "Spectra_Filepath", "Label", "Sample"])
        .to_csv(sep="\t", index=False, lineterminator="\n")
        + "\n"
        + pd.DataFrame(sample_rows, columns=["Sample", "MSstats_Condition", "MSstats_BioReplicate"])
        .to_csv(sep="\t", index=False, lineterminator="\n")
    )


def run_design_tsv(params: Mapping, design: Optional[pd.DataFrame] = None) -> str:
    """The design.tsv of a run, from its parameters (and its LFQ design table, if any)."""
    files = [str(f) for f in params.get("mzML-files", []) or []]
    if params.get("analysis-mode", "LFQ") == "LFQ":
        if design is None or design.empty:
            design = align_design(pd.DataFrame(columns=["file", "condition", "replicate", "fraction"]), files)
        else:
            design = align_design(design, files)
        return lfq_design_tsv(design, params)
    from src.common.isobaric import label_type
    return isobaric_design_tsv(files, params, label_type(dict(params)))


def write_design(results_dir, params: Mapping, design: Optional[pd.DataFrame] = None) -> Path:
    """Write the run's design.tsv into its combined/ folder and return the path."""
    path = combined_paths(results_dir)["design"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(run_design_tsv(params, design), encoding="utf-8")
    return path


def filtered_idxmls(results_dir) -> List[Path]:
    """The run's per-file idXMLs after the PSM FDR filter, in file order."""
    return sorted(Path(results_dir, "psm_filter").glob("*.idXML"))


def merge_command(idxmls: List, out) -> List[str]:
    """The IDMerger call that merges ``idxmls`` into ``out``, noting each PSM's source file."""
    return ["IDMerger", "-in", *[str(p) for p in idxmls], "-out", str(out), "-annotate_file_origin", "true"]


def build_missing(
    results_dir,
    params: Mapping,
    workflow_dir,
    run: Callable[[List[str]], bool] = None,
) -> Dict[str, bool]:
    """Create whichever of a finished run's combined files are missing.

    For runs made before these exports existed. The design comes from the
    design ProteomicsLFQ got (fractionated runs), else from the workflow's
    current design table aligned to the run's files; the conditions always
    come from the run's own parameters.
    """
    run = run or _run
    paths = combined_paths(results_dir)
    done = {}
    if not paths["idxml"].is_file():
        idxmls = filtered_idxmls(results_dir)
        if idxmls:
            paths["idxml"].parent.mkdir(parents=True, exist_ok=True)
            done["idxml"] = run(merge_command(idxmls, paths["idxml"])) and paths["idxml"].is_file()
    if not paths["design"].is_file():
        design = None
        for candidate in (Path(results_dir, "quant_results", "experimental-design.tsv"), Path(workflow_dir, "experimental-design.tsv")):
            if candidate.is_file():
                design = read_design(candidate)
                break
        write_design(results_dir, params, design)
        done["design"] = True
    return done


def load_run_params(params_file) -> dict:
    try:
        with open(params_file, "r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def _run(command: List[str]) -> bool:
    try:
        return subprocess.run(command, capture_output=True, check=False).returncode == 0
    except OSError:
        return False
