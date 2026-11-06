"""The identified MS2 spectra of a run, kept with its results.

A finished run writes, per input file, a small mzML with only the MS2 spectra
that have a PSM in its final identifications (``results/spectra/<file>.mzML``,
named like the input). The run's ZIP carries them, so an opened run shows its
PSM spectra when the original mzML files are not there.

Streamlit-free (pyOpenMS only), so it also runs in the RQ worker.
"""
from __future__ import annotations

import re
from pathlib import Path

import pyopenms as oms

# pyOpenMS >= 3.5.0 loads identifications into this container, older ones into a list
PeptideIdentificationList = getattr(oms, "PeptideIdentificationList", None)

SPECTRA_DIR = "spectra"
_SCAN = re.compile(r"scan=(\d+)")


def _scan(native_id: str) -> int:
    match = _SCAN.search(native_id)
    return int(match.group(1)) if match else 0


def identified_refs(idxml_path) -> tuple[set[str], set[int]]:
    """Spectrum references (native IDs) and scan numbers of the PSMs in an idXML.

    Only peptide identifications with at least one hit count.
    """
    proteins = []
    peptides = PeptideIdentificationList() if PeptideIdentificationList is not None else []
    oms.IdXMLFile().load(str(idxml_path), proteins, peptides)
    refs, scans = set(), set()
    for pep in peptides:
        if not pep.getHits() or not pep.metaValueExists("spectrum_reference"):
            continue
        ref = pep.getMetaValue("spectrum_reference")
        if isinstance(ref, bytes):
            ref = ref.decode()
        refs.add(str(ref))
        if _scan(str(ref)):
            scans.add(_scan(str(ref)))
    return refs, scans


def write_spectra_subset(mzml_path, idxml_paths, out_path) -> int:
    """Write the MS2 spectra of ``mzml_path`` identified in ``idxml_paths`` to ``out_path``.

    A spectrum is kept when its native ID is a PSM's spectrum reference, or
    its scan number is that of one (the viewers look spectra up by scan
    number). Spectra are copied unchanged, so native IDs, precursors and
    instrument settings stay as in the input; chromatograms are dropped.

    Returns:
        The number of spectra written.
    """
    refs, scans = set(), set()
    for idxml in idxml_paths:
        r, s = identified_refs(idxml)
        refs |= r
        scans |= s

    options = oms.PeakFileOptions()
    options.setMSLevels([2])
    exp = oms.MSExperiment()
    mzml = oms.MzMLFile()
    mzml.setOptions(options)
    mzml.load(str(mzml_path), exp)

    kept = [
        spec for spec in exp
        if spec.getMSLevel() == 2 and (spec.getNativeID() in refs or _scan(spec.getNativeID()) in scans)
    ]
    out = oms.MSExperiment(exp)  # keeps the run's instrument and sample settings
    out.setSpectra(kept)
    out.setChromatograms([])
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".part")
    writer = oms.MzMLFile()
    write_options = writer.getOptions()
    write_options.setCompression(True)  # zlib: lossless, keeps the run's ZIP small
    writer.setOptions(write_options)
    writer.store(str(tmp), out)
    tmp.replace(out_path)
    return len(kept)


def write_identified_spectra(mzml_files, idxml_files, results_dir) -> dict[str, int]:
    """One spectra subset per input file, under ``<results_dir>/spectra/``.

    ``idxml_files[i]`` holds the final identifications of ``mzml_files[i]``.
    Each subset has its input's file name, so the identifications still point
    at it.

    Returns:
        Spectra written per file name.
    """
    out_dir = Path(results_dir, SPECTRA_DIR)
    written = {}
    for mzml, idxml in zip(mzml_files, idxml_files):
        name = Path(mzml).name
        written[name] = write_spectra_subset(mzml, [idxml], out_dir / name)
    return written


def spectra_sources(mzml_files, results_dir) -> list[Path | None]:
    """Where to read each input's spectra: the input when it is there, else the run's subset.

    The order follows ``mzml_files``, which numbers the files in the viewer
    caches. ``None`` marks a file with neither.
    """
    sources = []
    for f in mzml_files:
        original = Path(f)
        subset = Path(results_dir, SPECTRA_DIR, original.name)
        if original.is_file():
            sources.append(original)
        elif subset.is_file():
            sources.append(subset)
        else:
            sources.append(None)
    return sources


def has_spectra(results_dir) -> bool:
    """Whether a run carries its identified spectra (runs made before this was added do not)."""
    return Path(results_dir, SPECTRA_DIR).is_dir()
