"""Why a run ended without peptide-spectrum matches, in words a user can act on.

The workflow stops when the search, the rescoring or the PSM FDR filter leaves
no identification. The results pages then have nothing to show; this module
reads what the run did write (idXML files per step, the input mzML files and
the run's parameters) and says which step emptied the result, at which FDR,
and what usually causes it. Streamlit-free, so it can be tested directly.
"""

from __future__ import annotations

from pathlib import Path

from src.common.fdr_display import fdr_applied, fdr_label

# Steps in run order: (results subfolder, what the step is called on screen).
# Runs from before ProSE wrote their search hits to comet_results.
LEGACY_FOLDERS = {"search_results": "comet_results"}
STEPS = (
    ("search_results", "database search"),
    ("percolator_results", "rescoring"),
    ("psm_filter", "PSM FDR filter"),
)
# ProteomicsLFQ's own default when neither the user nor the workflow set it
PROTEOMICSLFQ_PROTEIN_FDR = 0.05
# Below this many MS2 spectra the target-decoy estimate is coarse: with n target
# PSMs the smallest q-value one false match allows is about 1/n
FEW_SPECTRA = 1000
Q_VALUE = "MS:1001491"  # Percolator's PSM-level q-value
# What ProteomicsLFQ says when its protein FDR filter removes every protein
NO_PROTEINS = "No proteins left after FDR filtering"


def run_fdr_levels(params: dict) -> dict:
    """PSM and protein FDR of a run in percent; None where that filter did not run."""
    psm = min(float(params.get("psm-fdr-percent", 1.0)), 100.0)
    protein = None
    if params.get("analysis-mode", "LFQ") == "LFQ":
        value = params.get("ProteomicsLFQ", {}).get("proteinFDR")
        if value is None:
            value = params.get("_defaults", {}).get("ProteomicsLFQ", {}).get("proteinFDR")
        if value is None:
            value = PROTEOMICSLFQ_PROTEIN_FDR
        protein = min(float(value), 1.0) * 100.0
    # Isobaric runs infer proteins without a protein FDR filter
    return {
        "psm": psm if fdr_applied(psm) else None,
        "protein": protein if fdr_applied(protein) else None,
    }


def fdr_summary(params: dict) -> str:
    """'PSM FDR 1 % (q ≤ 0.01) · protein FDR not applied'."""
    levels = run_fdr_levels(params)
    return f"PSM FDR {fdr_label(levels['psm'])} · protein FDR {fdr_label(levels['protein'])}"


def _target_hits(idxml_path: Path) -> tuple[int, float | None]:
    """Target top hits of an idXML and the smallest q-value among them."""
    from src.common.results_helpers import load_idxml  # pyOpenMS only when needed

    n, best_q = 0, None
    for pep in load_idxml(str(idxml_path))[1]:
        hits = pep.getHits()
        if not hits:
            continue
        hit = hits[0]
        decoy = hit.getMetaValue("target_decoy") if hit.metaValueExists("target_decoy") else ""
        if isinstance(decoy, bytes):
            decoy = decoy.decode()
        if str(decoy).startswith("decoy"):
            continue
        n += 1
        if hit.metaValueExists(Q_VALUE):
            q = float(hit.getMetaValue(Q_VALUE))
            best_q = q if best_q is None else min(best_q, q)
    return n, best_q


def step_counts(results_dir: Path) -> list[dict]:
    """Target PSMs (and best q-value) after each step that wrote output."""
    counts = []
    for folder, label in STEPS:
        files = sorted(Path(results_dir, folder).glob("*.idXML"))
        if not files and folder in LEGACY_FOLDERS:
            files = sorted(Path(results_dir, LEGACY_FOLDERS[folder]).glob("*.idXML"))
        if not files:
            break
        n, best_q = 0, None
        for f in files:
            fn, fq = _target_hits(f)
            n += fn
            if fq is not None:
                best_q = fq if best_q is None else min(best_q, fq)
        counts.append({"step": label, "targets": n, "best_q": best_q})
    return counts


def _contains(path: Path, tag: bytes, chunk: int = 1 << 20) -> bool:
    """True when a file contains tag, read as bytes without parsing."""
    tail = b""
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                return False
            if tag in tail + block:
                return True
            tail = block[-len(tag):]


def filter_has_targets(results_dir: Path) -> bool:
    """Whether the PSM filter step kept a target PSM, mostly without parsing the idXML."""
    for f in Path(results_dir, "psm_filter").glob("*.idXML"):
        if _contains(f, b'name="target_decoy" value="target'):
            return True
        if _contains(f, b"<PeptideHit") and _target_hits(f)[0]:
            return True  # hits without a target/decoy label
    return False


def _log_mentions(run_dir: Path, text: str) -> bool:
    log = Path(run_dir, "logs", "minimal.log")
    return log.is_file() and _contains(log, text.encode())


def _no_proteins(results_dir: Path, params: dict) -> dict:
    """PSMs passed, but the protein FDR filter of ProteomicsLFQ kept no protein."""
    counts = step_counts(results_dir)
    levels = run_fdr_levels(params)
    passed = counts[-1]["targets"] if counts else 0
    counts.append({"step": "protein FDR filter", "targets": 0, "best_q": None, "unit": "proteins"})
    return {
        "step": "protein FDR filter",
        "what": "protein",
        "psm_fdr": fdr_label(levels["psm"]),
        "protein_fdr": fdr_label(levels["protein"]),
        "counts": counts,
        "ms2_spectra": None,
        "causes": [
            f"The {passed:,} PSMs that passed the PSM FDR do not add up to any protein at a "
            f"protein FDR of {fdr_label(levels['protein'])}: they are spread thinly over many "
            "proteins, as wrong matches are. The protein FDR is proteinFDR on the "
            "Quantification tab.",
            "The protein database does not fit the sample: wrong species, or a FASTA without "
            "the proteins you expect. Correct matches pile up on the proteins really present.",
            "Too few spectra: a small test file may identify a few peptides but not enough "
            "per protein to pass a protein FDR.",
        ],
    }


def diagnose_empty_run(run_dir: Path, params: dict, ms2_spectra=None) -> dict | None:
    """Explain a run whose PSM filter step left no target PSM; None when it has results.

    ms2_spectra is the number of MS2 spectra over the run's input files, when known,
    or a function returning it (only called for an empty run, since it reads the mzML).
    Returns the step that emptied the result, the FDR levels, the PSM count per
    step and a list of likely causes, most likely first.
    """
    results_dir = Path(run_dir, "results")
    if filter_has_targets(results_dir):
        if _log_mentions(run_dir, NO_PROTEINS):
            return _no_proteins(results_dir, params)
        return None
    counts = step_counts(results_dir)
    if not counts:
        return None  # never got to write a search result: not an empty result
    if counts[-1]["step"] == STEPS[-1][1] and counts[-1]["targets"] > 0:
        return None
    emptied = next((c for c in counts if c["targets"] == 0), None)
    if emptied is None:
        # Some step left PSMs but the run stopped before the filter: a crash, not this case
        return None

    if callable(ms2_spectra):
        ms2_spectra = ms2_spectra()
    levels = run_fdr_levels(params)
    psm = levels["psm"]
    causes = []
    few = ms2_spectra is not None and ms2_spectra < FEW_SPECTRA
    if emptied["step"] == "PSM FDR filter":
        before = next((c for c in counts if c["step"] == "rescoring"), counts[0])
        best_q = before["best_q"]
        detail = f"{before['targets']:,} target PSMs reached the filter"
        if best_q is not None:
            detail += f"; the best had q = {best_q:.2g}"
        causes.append(
            f"The PSM FDR of {fdr_label(psm)} is stricter than any match reached ({detail}). "
            "Try 1 % or 5 % on the False discovery rate tab."
        )
        if few:
            causes.append(
                f"Only {ms2_spectra:,} MS2 spectra: with few matches the target-decoy estimate "
                "cannot reach small q-values. Add more files or a longer run."
            )
        causes.append(
            "Matches that are mostly wrong also end here: check that the FASTA is the sample's "
            "species and that the mass tolerances and enzyme fit the instrument."
        )
    else:
        if few:
            causes.append(
                f"Only {ms2_spectra:,} MS2 spectra in the input: check that the files hold "
                "fragmentation (MS2) data and were not filtered or cut short."
            )
        causes.append(
            "The protein database does not fit the sample: wrong species, or a FASTA without "
            "the proteins you expect."
        )
        causes.append(
            "The precursor or fragment mass tolerance, or the enzyme, does not fit the data "
            "(e.g. a high-resolution setting on ion-trap spectra). Pick your Instrument in "
            "Your analysis on the Describe experiment page."
        )
        if not few:
            causes.append(
                "Too few spectra to estimate an FDR: very small test files often end here."
            )
    return {
        "step": emptied["step"],
        "what": "PSM",
        "psm_fdr": fdr_label(psm),
        "protein_fdr": fdr_label(levels["protein"]),
        "counts": counts,
        "ms2_spectra": ms2_spectra,
        "causes": causes,
    }
