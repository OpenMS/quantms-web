"""Rescore the matches of all files together and split the result per file.

Percolator learns the most from a study when it sees every file's matches at
once: one model, and q-values and posterior error probabilities that mean the
same in every file. PercolatorAdapter does this when it is given all search
results as one ``-in`` list, but it then writes a single idXML. The steps after
it (FDR filtering, quantification, the result views) work file by file, so
``split_pooled`` writes each file's matches back into an idXML of its own.

PercolatorAdapter (OpenMS 3.6) marks where each match came from: with more than
one input, the spectrum reference of every PeptideIdentification becomes
``file=<i>,<reference>``, ``<i>`` being the position of its file in the
``-in`` list. Its protein runs are merged into one and every match points at
it. ``split_pooled`` groups the matches by that position, restores their
spectrum reference, and gives each file its own protein run (the one of its
search result, with the annotations Percolator added to the merged run).

Streamlit-free, so it runs in the workflow process.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import List, Sequence

import pyopenms as oms

_FILE_PREFIX = re.compile(r"^file=(\d+),(.*)$")

PeptideIdentificationList = getattr(oms, "PeptideIdentificationList", None)


def _peptide_list():
    return PeptideIdentificationList() if PeptideIdentificationList is not None else []


def _load(path: str):
    proteins: list = []
    peptides = _peptide_list()
    oms.IdXMLFile().load(str(path), proteins, peptides)
    return proteins, peptides


def _original_reference(reference: str, pep) -> str:
    """The spectrum reference a match had before Percolator added its file.

    Percolator falls back to ``scan=<spectrum_id>`` or ``index=<n>`` for a match
    without a spectrum reference; such a match had none.
    """
    if reference.startswith("index="):
        return ""
    if reference.startswith("scan=") and pep.metaValueExists("spectrum_id"):
        if reference == f"scan={pep.getMetaValue('spectrum_id')}":
            return ""
    return reference


def split_pooled(merged: str, inputs: Sequence[str], outputs: Sequence[str]) -> List[int]:
    """Write the matches of ``merged`` from input ``i`` to ``outputs[i]``.

    Returns the number of PeptideIdentifications written per file. Raises
    ``ValueError`` when ``inputs`` and ``outputs`` differ in length or a match
    carries no file position (a result Percolator was not given a list for).
    """
    if len(inputs) != len(outputs):
        raise ValueError(f"{len(inputs)} inputs but {len(outputs)} outputs")
    merged_proteins, merged_peptides = _load(merged)
    merged_run = merged_proteins[0] if merged_proteins else None

    groups: List[list] = [[] for _ in inputs]
    for pep in merged_peptides:
        reference = pep.getSpectrumReference()
        if isinstance(reference, bytes):
            reference = reference.decode()
        match = _FILE_PREFIX.match(reference)
        if match is None:
            raise ValueError(
                f"Match without a file position in {Path(merged).name}: {reference!r}"
            )
        index = int(match.group(1))
        if index >= len(inputs):
            raise ValueError(f"Match from file {index}, but only {len(inputs)} files were rescored")
        pep.setSpectrumReference(_original_reference(match.group(2), pep))
        groups[index].append(pep)

    counts = []
    for source, target, peptides in zip(inputs, outputs, groups):
        proteins, _ = _load(source)
        run = proteins[0]
        if merged_run is not None:
            keys: list = []
            merged_run.getKeys(keys)
            for key in keys:
                if not run.metaValueExists(key):
                    run.setMetaValue(key, merged_run.getMetaValue(key))
        identifier = run.getIdentifier()
        container = _peptide_list()
        for pep in peptides:
            pep.setIdentifier(identifier)
            if PeptideIdentificationList is not None:
                container.push_back(pep)
            else:
                container.append(pep)
        Path(target).parent.mkdir(parents=True, exist_ok=True)
        oms.IdXMLFile().store(str(target), [run], container)
        counts.append(len(peptides))
    return counts
