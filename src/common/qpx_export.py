"""QPX (quantms Parquet) export of a run: the PSM view and the provenance view.

QPX (https://github.com/bigbio/qpx, the successor of quantms.io) stores a
proteomics result as Parquet views. This module writes two of them for a run:

- ``qpx/<run>.psm.parquet``: every PSM after the PSM FDR filter, one row per
  spectrum (its top hit), following QPX's ``psm`` schema (spec 1.2).
- ``qpx/<run>.provenance.parquet``: one row per processing step, with the
  tool, its version, the parameters the workflow set on its command line and
  the tool's full configuration (every parameter) as JSON.

OpenMS 3.5 has no writer for this format (QuantmsIOConverter writes the older
quantms.io 1.0 layout, without protein accessions and with the native ID in
place of the run name), so the PSMs are read with pyOpenMS and written with
pyarrow. The schemas below are transcribed field by field from QPX's
``psm.yaml``, ``provenance.yaml`` and ``types.yaml``.

The provenance is read back from what the run did: the commands it logged
(every effective parameter is an explicit command-line argument) and each
tool's ``-write_ini`` (all its parameters and the OpenMS version). This module
makes no Streamlit calls, so it runs in a queue worker too.
"""
from __future__ import annotations

import hashlib
import json
import re
import shlex
import subprocess
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List, Mapping, Optional

import pyarrow as pa
import pyarrow.parquet as pq

QPX_SPEC_VERSION = "1.2"
QPX_DIR = "qpx"
PSM_IDENTITY = ("peptidoform", "charge", "run_file_name", "scan")

# Shared QPX struct types (types.yaml)
SCORE = pa.struct([
    pa.field("score_name", pa.string()),
    pa.field("score_value", pa.float64()),
    pa.field("higher_better", pa.bool_(), nullable=True),
])
CV_PARAM = pa.struct([pa.field("cv_name", pa.string()), pa.field("cv_value", pa.string())])
MODIFICATION_POSITION = pa.struct([
    pa.field("position", pa.int32()),
    pa.field("amino_acid", pa.string(), nullable=True),
    pa.field("scores", pa.list_(SCORE), nullable=True),
])
MODIFICATION = pa.struct([
    pa.field("name", pa.string()),
    pa.field("accession", pa.string(), nullable=True),
    pa.field("positions", pa.list_(MODIFICATION_POSITION)),
])
PARAMETER = pa.struct([pa.field("key", pa.string()), pa.field("value", pa.string())])

# psm.yaml: required fields are not nullable unless the spec says so
PSM_SCHEMA = pa.schema([
    pa.field("psm_id", pa.int64(), nullable=False),
    pa.field("sequence", pa.string(), nullable=False),
    pa.field("peptidoform", pa.string(), nullable=False),
    pa.field("modifications", pa.list_(MODIFICATION)),
    pa.field("charge", pa.int16(), nullable=False),
    pa.field("posterior_error_probability", pa.float64()),
    pa.field("is_decoy", pa.bool_(), nullable=False),
    pa.field("calculated_mz", pa.float32()),
    pa.field("observed_mz", pa.float32(), nullable=False),
    pa.field("mass_error_ppm", pa.float32()),
    pa.field("additional_scores", pa.list_(SCORE)),
    pa.field("predicted_rt", pa.float32()),
    pa.field("run_file_name", pa.string(), nullable=False),
    pa.field("cv_params", pa.list_(CV_PARAM)),
    pa.field("scan", pa.list_(pa.int32()), nullable=False),
    pa.field("rt", pa.float32()),
    pa.field("missed_cleavages", pa.int16()),
    pa.field("protein_accessions", pa.list_(pa.string())),
])

# provenance.yaml
PROVENANCE_SCHEMA = pa.schema([
    pa.field("step_order", pa.int32(), nullable=False),
    pa.field("step_category", pa.string(), nullable=False),
    pa.field("step_name", pa.string(), nullable=False),
    pa.field("tool_name", pa.string(), nullable=False),
    pa.field("tool_version", pa.string()),
    pa.field("tool_uri", pa.string()),
    pa.field("parameters", pa.list_(PARAMETER)),
    pa.field("config", pa.string()),
    pa.field("output_views", pa.list_(pa.string())),
])

# What each TOPP tool of the workflow does, as (step_category, step_name)
STEPS = {
    "DecoyDatabase": ("database", "decoy_database"),
    "CometAdapter": ("database_search", "sequence_search"),
    "PercolatorAdapter": ("rescoring", "psm_rescoring"),
    "IDFilter": ("filtering", "fdr_filtering"),
    "IDMerger": ("format_conversion", "identification_merging"),
    "ProteomicsLFQ": ("quantification", "label_free_quantification"),
    "IsobaricAnalyzer": ("quantification", "isobaric_quantification"),
    "IDMapper": ("quantification", "identification_mapping"),
    "FileMerger": ("quantification", "map_merging"),
    "ProteinInference": ("protein_inference", "protein_inference"),
    "IDConflictResolver": ("protein_inference", "conflict_resolution"),
    "ProteinQuantifier": ("quantification", "protein_quantification"),
}
# The search engines behind adapters: (executable parameter, default, version pattern)
ENGINES = {
    "CometAdapter": ("comet_executable", "comet", r'Comet version "([^"]+)"', ["comet"]),
    "PercolatorAdapter": ("percolator_executable", "percolator", r"Percolator version ([\d.]+\d)", ["percolator", "-h"]),
}
TOPP_URI = "https://openms.de/doxygen/release/latest/html/TOPP_{tool}.html"


# ----------------------------------------------------------------------------- PSMs

def derive_id(values: list) -> int:
    """QPX's opaque identity: BLAKE2b-64 of the canonical JSON of the composite, signed."""
    encoded = json.dumps(list(values), separators=(",", ":"), sort_keys=True, ensure_ascii=True).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(encoded, digest_size=8).digest(), "big", signed=True)


def scan_number(native_id: str) -> Optional[int]:
    """The scan number in a native ID (``... scan=1234``, ``index=12`` or a bare number)."""
    m = re.search(r"(?:^|\s)scan=(\d+)", native_id or "") or re.search(r"(?:^|\s)index=(\d+)", native_id or "")
    if m:
        return int(m.group(1))
    return int(native_id) if (native_id or "").isdigit() else None


def missed_cleavages(sequence: str) -> int:
    """Internal K/R not followed by P (trypsin's rule)."""
    return sum(1 for i, aa in enumerate(sequence[:-1]) if aa in "KR" and sequence[i + 1] != "P")


def modifications(aa_sequence) -> Optional[list]:
    """QPX modifications of a pyOpenMS AASequence: name, UNIMOD accession, 1-based positions.

    N-terminal modifications are at position 0 and C-terminal ones at length + 1.
    """
    found = {}

    def add(mod, position, residue):
        name = mod.getId() or mod.getFullId()
        unimod = re.search(r"(\d+)$", mod.getUniModAccession() or "")
        accession = f"UNIMOD:{unimod.group(1)}" if unimod else None
        entry = found.setdefault(name, {"name": name, "accession": accession, "positions": []})
        entry["positions"].append({"position": position, "amino_acid": residue, "scores": None})

    if aa_sequence.hasNTerminalModification():
        add(aa_sequence.getNTerminalModification(), 0, "N-term")
    for i in range(aa_sequence.size()):
        residue = aa_sequence.getResidue(i)
        if residue.isModified():
            add(residue.getModification(), i + 1, residue.getOneLetterCode())
    if aa_sequence.hasCTerminalModification():
        add(aa_sequence.getCTerminalModification(), aa_sequence.size() + 1, "C-term")
    return list(found.values()) or None


def _meta(hit, key):
    if not hit.metaValueExists(key):
        return None
    value = hit.getMetaValue(key)
    return value.decode() if isinstance(value, bytes) else value


def _float(value) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# Scores kept in additional_scores: (meta value or "main", QPX score name, higher is better)
SCORES = [
    ("MS:1001491", "percolator q-value", False),
    ("MS:1001492", "percolator score", True),
    ("COMET:xcorr", "comet xcorr", True),
    ("COMET:deltaCn", "comet deltacn", True),
    ("expect", "comet expectation value", False),
    ("q-value", "q-value", False),
]


def psm_rows(idxml_path, run_file_name: str, to_proforma) -> List[dict]:
    """The QPX psm rows of one idXML (top hit of each spectrum)."""
    from src.common.results_helpers import load_idxml

    _, peptides = load_idxml(idxml_path)
    rows = []
    for pep in peptides:
        hits = pep.getHits()
        if not hits:
            continue
        hit = hits[0]
        seq = hit.getSequence()
        charge = int(hit.getCharge())
        peptidoform = to_proforma(seq.toString())
        scan = scan_number(str(pep.getMetaValue("spectrum_reference")) if pep.metaValueExists("spectrum_reference") else "")
        scans = [scan] if scan is not None else []
        calculated = seq.getMZ(charge) if charge else None
        observed = pep.getMZ()
        score_type = pep.getScoreType()
        pep_value = _float(_meta(hit, "MS:1001493"))
        if pep_value is None and "Posterior Error Probability" in score_type:
            pep_value = hit.getScore()
        scores = [
            {"score_name": name, "score_value": float(v), "higher_better": higher}
            for key, name, higher in SCORES
            if (v := _float(_meta(hit, key))) is not None
        ]
        if "Posterior Error Probability" not in score_type and not any(s["score_name"] == score_type for s in scores):
            scores.insert(0, {"score_name": score_type, "score_value": float(hit.getScore()),
                              "higher_better": bool(pep.isHigherScoreBetter())})
        unmodified = seq.toUnmodifiedString()
        rows.append({
            "psm_id": derive_id([peptidoform, charge, run_file_name, scans]),
            "sequence": unmodified,
            "peptidoform": peptidoform,
            "modifications": modifications(seq),
            "charge": charge,
            "posterior_error_probability": pep_value,
            "is_decoy": str(_meta(hit, "target_decoy") or "").startswith("decoy"),
            "calculated_mz": calculated,
            "observed_mz": observed,
            "mass_error_ppm": 1e6 * (observed - calculated) / calculated if calculated else None,
            "additional_scores": scores or None,
            "predicted_rt": None,
            "run_file_name": run_file_name,
            "cv_params": None,
            "scan": scans,
            "rt": pep.getRT(),
            "missed_cleavages": missed_cleavages(unmodified),
            "protein_accessions": sorted({ev.getProteinAccession() for ev in hit.getPeptideEvidences()}) or None,
        })
    return rows


def footer(file_type: str, extra: Optional[Mapping[str, str]] = None) -> dict:
    meta = {
        "qpx_version": QPX_SPEC_VERSION,
        "file_type": file_type,
        "creator": "OpenDDA",
        "creation_date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "uuid": str(uuid.uuid4()),
        "compression_format": "zstd",
    }
    meta.update(extra or {})
    return meta


def write_table(rows: List[dict], schema: pa.Schema, path, metadata: Mapping[str, str]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pylist(rows, schema=schema).replace_schema_metadata(dict(metadata))
    pq.write_table(table, str(path), compression="zstd")
    return path


def write_psms(idxmls: Iterable, mzml_names: List[str], path, to_proforma) -> int:
    """Write the PSM view of ``idxmls`` and return its number of rows."""
    from src.common.results_helpers import extract_filename_from_idxml

    rows = []
    for idxml in idxmls:
        mzml = extract_filename_from_idxml(Path(idxml), mzml_names)
        run = mzml[:-5] if mzml.lower().endswith(".mzml") else Path(mzml).stem
        rows.extend(psm_rows(idxml, run, to_proforma))
    write_table(rows, PSM_SCHEMA, path, footer("psm_file", {
        "scan_format": "scan",
        "primary_key": "psm_id",
        "identity_composite": ",".join(PSM_IDENTITY),
    }))
    return len(rows)


# ----------------------------------------------------------------------------- provenance

def logged_commands(log_text: str) -> List[List[str]]:
    """The commands a run finished, in order, from its commands-and-run-times log."""
    commands = []
    lines = log_text.splitlines()
    for i, line in enumerate(lines):
        if line.strip() == "Process finished:" and i + 1 < len(lines):
            try:
                commands.append(shlex.split(lines[i + 1]))
            except ValueError:
                commands.append(lines[i + 1].split())
    return commands


def parse_command(command: List[str], known: Iterable[str]) -> dict:
    """``{parameter: value}`` of a TOPP command line; ``known`` are the tool's parameter names.

    A token is a parameter name only when the tool has a parameter of that
    name, so negative numbers stay values. Several values are joined by a space.
    """
    known = set(known)
    params, key = {}, None
    for token in command[1:]:
        if token.startswith("-") and token[1:] in known:
            key = token[1:]
            params[key] = []  # given twice: the tool takes the last one
        elif key is not None:
            params[key].append(token)
    return {k: " ".join(v) if v else "true" for k, v in params.items()}


def tool_ini(tool: str, run=subprocess.run) -> tuple:
    """(OpenMS version, {parameter: default}) of a TOPP tool, from its ``-write_ini``."""
    import pyopenms as poms

    with tempfile.TemporaryDirectory() as tmp:
        ini = Path(tmp, f"{tool}.ini")
        try:
            run([tool, "-write_ini", str(ini)], capture_output=True, check=False, timeout=60)
        except (OSError, subprocess.SubprocessError):
            return None, {}
        if not ini.is_file():
            return None, {}
        text = ini.read_text(errors="ignore")
        version = re.search(r'name="version" value="([^"]+)"', text)
        param = poms.Param()
        poms.ParamXMLFile().load(str(ini), param)
    defaults = {}
    prefix = f"{tool}:1:"
    for key in param.keys():
        name = key.decode() if isinstance(key, bytes) else str(key)
        if not name.startswith(prefix):
            continue
        value = param.getValue(key)
        if isinstance(value, bytes):
            value = value.decode()
        elif isinstance(value, list):
            value = " ".join(v.decode() if isinstance(v, bytes) else str(v) for v in value)
        defaults[name[len(prefix):]] = value
    return (version.group(1) if version else None), defaults


def engine_version(tool: str, params: Mapping[str, str], run=subprocess.run) -> Optional[str]:
    if tool not in ENGINES:
        return None
    key, default, pattern, args = ENGINES[tool]
    exe = params.get(key) or default
    try:
        out = run([exe, *args[1:]], capture_output=True, text=True, check=False, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    m = re.search(pattern, (out.stdout or "") + (out.stderr or ""))
    return m.group(1) if m else None


def _relative(value: str, workflow_dir) -> str:
    """Paths inside the workflow folder, relative to it (no server paths in the export)."""
    root = str(Path(workflow_dir))
    parts = []
    for item in str(value).split(" "):
        if root and root in item:
            item = item.split(root, 1)[1].lstrip("/\\")
        parts.append(item)
    return " ".join(parts)


def provenance_rows(log_text: str, workflow_dir, params: Mapping, app_version: str,
                    ini=tool_ini, engine=engine_version) -> List[dict]:
    """One provenance row per tool the run used (first call of each), then the export itself."""
    rows, by_tool = [], {}
    for command in logged_commands(log_text):
        tool = Path(command[0]).name if command else ""
        if tool:
            by_tool.setdefault(tool, []).append(command)
    for tool, commands in by_tool.items():
        version, defaults = ini(tool)
        if not defaults and tool not in STEPS:
            continue  # not a TOPP tool (e.g. a Python step)
        # One row per tool; a value that differs between its calls (each file's
        # input and output) lists every call's value
        calls = [parse_command(c, defaults) for c in commands]
        given = {}
        for call in calls:
            for k in call:
                values = list(dict.fromkeys(_relative(c[k], workflow_dir) for c in calls if k in c))
                given[k] = "; ".join(values)
        config = {k: str(v) for k, v in defaults.items()}
        config.update(given)
        tool_version = version
        engine_v = engine(tool, given)
        if engine_v:
            given["engine_version"] = engine_v
        category, step = STEPS.get(tool, ("processing", tool.lower()))
        rows.append({
            "step_category": category,
            "step_name": step,
            "tool_name": f"OpenMS/{tool}",
            "tool_version": tool_version,
            "tool_uri": TOPP_URI.format(tool=tool),
            "parameters": [{"key": k, "value": str(v)} for k, v in given.items()],
            "config": json.dumps(config, sort_keys=True),
            "output_views": ["psm"] if tool == "IDFilter" else [],
        })
    database = next(
        (p["value"] for r in rows if r["tool_name"] == "OpenMS/CometAdapter" for p in r["parameters"] if p["key"] == "database"),
        None,
    )
    workflow = {
        "analysis_mode": params.get("analysis-mode", "LFQ"),
        "fasta_file": _relative(params.get("fasta-file", ""), workflow_dir),
        "search_database": database or "",
        "psm_fdr_percent": params.get("psm-fdr-percent", ""),
        "decoys_generated": params.get("generate-decoys", ""),
        "contaminants_added": params.get("add-contaminants", ""),
        "mzML_files": " ".join(Path(str(f)).name for f in params.get("mzML-files", []) or []),
    }
    rows.insert(0, {
        "step_category": "workflow",
        "step_name": "workflow",
        "tool_name": "OpenDDA",
        "tool_version": app_version,
        "tool_uri": "https://github.com/OpenMS/quantms-web",
        "parameters": [{"key": k, "value": str(v)} for k, v in workflow.items()],
        "config": json.dumps({k: v for k, v in params.items() if not k.startswith("_")}, sort_keys=True, default=str),
        "output_views": [],
    })
    rows.append({
        "step_category": "format_conversion",
        "step_name": "qpx_export",
        "tool_name": "OpenDDA",
        "tool_version": app_version,
        "tool_uri": "https://github.com/OpenMS/quantms-web",
        "parameters": [{"key": "qpx_version", "value": QPX_SPEC_VERSION}],
        "config": None,
        "output_views": ["psm", "provenance"],
    })
    for i, row in enumerate(rows, start=1):
        row["step_order"] = i
    return rows


def provenance_table(rows: List[dict]):
    """The provenance as a flat table for display: one row per step and parameter."""
    import pandas as pd

    flat = []
    for row in rows:
        config = json.loads(row["config"]) if row.get("config") and row["tool_name"] != "OpenDDA" else {}
        given = {p["key"]: p["value"] for p in row.get("parameters") or []}
        keys = list(given) + [k for k in config if k not in given]
        for key in keys or [""]:
            flat.append({
                "step": row["step_order"],
                "tool": row["tool_name"],
                "version": row.get("tool_version") or "",
                "parameter": key,
                "value": str(given.get(key, config.get(key, ""))),
                "set by": "workflow" if key in given else "tool default",
            })
    return pd.DataFrame(flat, columns=["step", "tool", "version", "parameter", "value", "set by"])


# ----------------------------------------------------------------------------- export

def qpx_paths(results_dir, name: str = "run") -> dict:
    folder = Path(results_dir, QPX_DIR)
    return {"psm": folder / f"{name}.psm.parquet", "provenance": folder / f"{name}.provenance.parquet"}


def export_run(results_dir, workflow_dir, params: Mapping, log_text: str, app_version: str,
               to_proforma=None, name: str = "run") -> dict:
    """Write a run's QPX psm and provenance views. Returns the paths and the PSM count."""
    if to_proforma is None:
        from src.common.peptide_format import to_proforma
    paths = qpx_paths(results_dir, name)
    idxmls = sorted(Path(results_dir, "psm_filter").glob("*.idXML"))
    mzml_names = [Path(str(f)).name for f in params.get("mzML-files", []) or []]
    n = write_psms(idxmls, mzml_names, paths["psm"], to_proforma)
    write_table(provenance_rows(log_text, workflow_dir, params, app_version), PROVENANCE_SCHEMA,
                paths["provenance"], footer("provenance_file"))
    return {**paths, "psms": n}


def app_version() -> str:
    """The app's version from settings.json (empty when it cannot be read)."""
    try:
        with open("settings.json", "r", encoding="utf-8") as f:
            return str(json.load(f).get("version", ""))
    except (OSError, ValueError):
        return ""


def read_provenance(path) -> List[dict]:
    return pq.read_table(str(path)).to_pylist()
