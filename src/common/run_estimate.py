"""A rough time and memory estimate for a run, made before it starts.

The numbers come from the inputs alone, read cheaply: the spectrum count an
indexed mzML states in its header (else its file size), the FASTA's size and
protein count, and the search settings. The constants were measured with
OpenMS 3.5 and Comet 2026.01 on the demo data (1,810 spectra per file, a
12,000-protein target+decoy database, 2 variable modifications):

- Comet holds every spectrum of a batch in memory, binned at the fragment
  tolerance: about 0.28 MB per MS2 spectrum at high resolution, on top of
  about 60 MB (183 MB with spectrum_batch_size 500, 550 MB for the whole file;
  the database size barely matters). The workflow searches up to max_threads
  files at once, each in its own Comet process.
- Comet needs about 0.032 CPU-seconds per spectrum against that database,
  roughly in proportion to the database size to the power 0.8 (4x the
  database took 3.1x as long); semi-specific digestion, more variable
  modifications and wide precursor windows multiply it.
- ProteomicsLFQ loads one mzML at a time; its peak is about 100 MB plus four
  times the largest mzML file.

An estimate is meant to tell a 5-minute run from a 5-hour one and a run that
fits in memory from one that does not, not to be exact.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Mapping, Optional

MB = 1024 ** 2
GB = 1024 ** 3

# Measured constants (see the module docstring)
COMET_BASE_MB = 60
COMET_MB_PER_SPECTRUM = 0.28          # high resolution (fragment bins of 0.02 Da)
COMET_MB_PER_SPECTRUM_LOW_RES = 0.03  # low resolution (bins of 1 Da)
COMET_CPU_S_PER_SPECTRUM = 0.032      # per MS2 spectrum against the reference database
REFERENCE_DB_BYTES = 7.4 * MB         # the target+decoy demo database those numbers used
LFQ_BASE_MB = 100
LFQ_MZML_FACTOR = 4
MS2_FRACTION = 0.9                    # share of MS2 spectra in a typical DDA file
SPECTRA_PER_MB = 120                  # when the mzML header states no count (demo: 1810 / 15 MB)
PER_FILE_OVERHEAD_S = 35              # Percolator, filtering and caches per file (15-85 s measured)
LFQ_S_PER_FILE = 5


@dataclass
class Estimate:
    files: int
    spectra: int
    proteins: int
    parallel: int
    minutes: float
    peak_gb: float
    comet_gb_per_file: float
    lfq_gb: float
    total_gb: Optional[float]
    search_factor: float
    advice: List[str] = field(default_factory=list)
    level: str = "ok"  # "ok", "long" or "high"


def spectrum_count(mzml: Path) -> int:
    """Spectra in an mzML: the count its spectrumList states near the top, else a guess from its size."""
    try:
        with open(mzml, "rb") as f:
            head = f.read(1 << 20)
    except OSError:
        return 0
    m = re.search(rb'<spectrumList[^>]*\scount="(\d+)"', head)
    if m:
        return int(m.group(1))
    return int(Path(mzml).stat().st_size / MB * SPECTRA_PER_MB)


def fasta_stats(fasta: Path) -> tuple:
    """(proteins, bytes) of a FASTA file."""
    proteins = 0
    try:
        with open(fasta, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 22), b""):
                proteins += chunk.count(b">")
        return proteins, Path(fasta).stat().st_size
    except OSError:
        return 0, 0


def _count_mods(value) -> int:
    if isinstance(value, list):
        return len([v for v in value if str(v).strip()])
    return len([v for v in str(value or "").split("\n") if v.strip()])


def search_factor(comet: Mapping) -> float:
    """How much longer than the reference settings a Comet search takes."""
    factor = 1.0
    termini = str(comet.get("num_enzyme_termini", "fully"))
    factor *= {"fully": 1.0, "semi": 8.0, "C-term unspecific": 8.0, "N-term unspecific": 8.0}.get(termini, 1.0)
    if str(comet.get("enzyme", "Trypsin")).lower() in ("unspecific cleavage", "no cleavage"):
        factor *= 50.0
    variable = _count_mods(comet.get("variable_modifications", "Oxidation (M)\nAcetyl (Protein N-term)"))
    max_mods = int(float(comet.get("max_variable_mods_in_peptide", 3) or 3))
    # each variable modification beyond the reference two adds candidates; more mods per peptide too
    factor *= 1.6 ** max(0, variable - 2) * 1.3 ** max(0, min(max_mods, variable) - 2)
    missed = int(float(comet.get("missed_cleavages", 2) or 2))
    factor *= 1.4 ** (missed - 2)
    tol = float(comet.get("precursor_mass_tolerance", 20.0) or 20.0)
    unit = str(comet.get("precursor_error_units", "ppm"))
    tol_da = tol if unit == "Da" else tol * 1e-6 * 1000  # at about 1,000 Da
    if tol_da > 0.02:
        factor *= min(400.0, (tol_da / 0.02) ** 0.6)
    return factor


def estimate(
    mzml_files: List[Path],
    fasta: Optional[Path],
    comet: Mapping,
    max_threads: int,
    decoys: bool = True,
    total_memory: Optional[int] = None,
) -> Estimate:
    """Time and memory a run would need, with advice when either is high."""
    sizes = [Path(f).stat().st_size if Path(f).is_file() else 0 for f in mzml_files]
    spectra = [spectrum_count(Path(f)) for f in mzml_files]
    ms2 = [s * MS2_FRACTION for s in spectra]
    proteins, db_bytes = fasta_stats(Path(fasta)) if fasta else (0, 0)
    if decoys:
        db_bytes *= 2
    n = len(mzml_files)
    parallel = max(1, min(n, max_threads))
    threads_each = max(1, max_threads // parallel)

    batch = int(float(comet.get("spectrum_batch_size", 20000) or 0))
    low_res = str(comet.get("instrument", "high_res")) == "low_res"
    per_spectrum = COMET_MB_PER_SPECTRUM_LOW_RES if low_res else COMET_MB_PER_SPECTRUM
    in_batch = [min(s, batch) if batch > 0 else s for s in ms2]
    largest_batch = max(in_batch, default=0)
    comet_mb = COMET_BASE_MB + per_spectrum * largest_batch
    lfq_mb = LFQ_BASE_MB + LFQ_MZML_FACTOR * max(sizes, default=0) / MB
    peak_mb = max(parallel * comet_mb, lfq_mb)

    factor = search_factor(comet)
    db_scale = (db_bytes / REFERENCE_DB_BYTES) ** 0.8 if db_bytes else 1.0
    cpu_s = sum(ms2) * COMET_CPU_S_PER_SPECTRUM * db_scale * factor
    # files are searched max_threads at a time, threads split between them
    search_s = cpu_s / max(1, parallel * threads_each)
    seconds = search_s + PER_FILE_OVERHEAD_S * n / parallel + LFQ_S_PER_FILE * n + 20

    est = Estimate(
        files=n, spectra=int(sum(spectra)), proteins=proteins, parallel=parallel,
        minutes=seconds / 60, peak_gb=peak_mb / 1024, comet_gb_per_file=comet_mb / 1024,
        lfq_gb=lfq_mb / 1024, total_gb=(total_memory / GB) if total_memory else None,
        search_factor=factor,
    )
    est.advice, est.level = advise(est, comet, batch, max_threads)
    return est


def advise(est: Estimate, comet: Mapping, batch: int, max_threads: int) -> tuple:
    """Advice for a run that needs much memory or time, and how serious it is."""
    advice, level = [], "ok"
    if est.total_gb and est.peak_gb > 0.7 * est.total_gb:
        level = "high"
        comet_part = est.parallel * est.comet_gb_per_file
        if comet_part >= est.lfq_gb:
            budget_mb = 0.6 * est.total_gb * 1024 / est.parallel
            suggested = int(max(500, (budget_mb - COMET_BASE_MB) / COMET_MB_PER_SPECTRUM) // 500 * 500)
            advice.append(
                f"Lower Comet's **spectrum_batch_size** (now {batch or 'all spectra'}) to about {suggested:,} "
                "on **Describe experiment**, tab *Database search*, with *Show advanced parameters* on. "
                "Comet then holds fewer spectra in memory at once; the results do not change."
            )
            if est.parallel > 1:
                advice.append(
                    f"Use fewer threads (**Settings** in the sidebar, now {max_threads}): {est.parallel} files "
                    "are searched at the same time and each needs its own memory."
                )
        else:
            advice.append(
                "The quantification step loads one whole mzML file at a time. Split very large "
                "files, or run the analysis on a machine with more memory."
            )
    if est.minutes > 120:
        level = "high" if level == "high" else "long"
        if est.files > 4:
            advice.append(
                "Try the settings on **two or three files** first: a short run shows whether "
                "they suit the data before you search all files."
            )
        if est.search_factor > 4:
            advice.append(
                "The search settings make the search much slower than a standard tryptic search "
                f"(about {est.search_factor:.0f}x): check semi-specific digestion, the number of variable "
                "modifications, missed cleavages and a wide precursor window."
            )
        advice.append(
            "Search only the proteome of your organism (plus contaminants): a larger database "
            "takes longer and costs identifications at the same FDR."
        )
    return advice, level


def total_memory() -> Optional[int]:
    """Bytes of RAM on this machine (None if unknown)."""
    try:
        import psutil
        return int(psutil.virtual_memory().total)
    except Exception:  # psutil missing or unsupported platform
        try:
            return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
        except (ValueError, OSError, AttributeError):
            return None


def format_minutes(minutes: float) -> str:
    if minutes < 1.5:
        return "about 1 minute"
    if minutes < 90:
        return f"about {round(minutes / 5) * 5 if minutes > 15 else round(minutes)} minutes"
    hours = minutes / 60
    return f"about {hours:.1f} hours" if hours < 10 else f"about {round(hours)} hours"
