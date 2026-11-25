"""QPX (quantms Parquet) export of a run: its PSM, protein group, run, sample and provenance views.

QPX (https://github.com/bigbio/qpx, the successor of quantms.io) stores a
proteomics result as Parquet views. This module writes these for a run:

- ``qpx/<run>.psm.parquet``: every PSM after the PSM FDR filter, one row per
  spectrum (its top hit), following QPX's ``psm`` schema (spec 1.2).
- ``qpx/<run>.provenance.parquet``: one row per processing step, with the
  tool, its version, the parameters the workflow set on its command line and
  the tool's full configuration (every parameter) as JSON.
- ``qpx/<run>.sample.parquet``: one row per sample (an isobaric run's
  channels are samples), with its condition and replicate and, from the SDRF
  imported with the data, its organism and other characteristics. Without an
  SDRF the organism is the one most entries of the protein database name.

OpenMS 3.5 has no writer for this format (QuantmsIOConverter writes the older
quantms.io 1.0 layout, without protein accessions and with the native ID in
place of the run name), so the PSMs are read with pyOpenMS and written with
pyarrow. The schemas below are transcribed field by field from QPX's
``psm.yaml``, ``provenance.yaml``, ``sample.yaml`` and ``types.yaml``.

The provenance is read back from what the run did: the commands it logged
(every effective parameter is an explicit command-line argument) and each
tool's ``-write_ini`` (all its parameters and the OpenMS version). This module
makes no Streamlit calls, so it runs in a queue worker too.
"""
from __future__ import annotations

import csv
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
# Mass difference between the 13C and 12C isotopes: a search engine that matches
# a precursor picked one or more isotope peaks off reports the offset per hit
C13_DELTA = 1.0033548378

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

# pg.yaml (the fields OpenDDA fills; the others are optional)
PEPTIDE_PER_PROTEIN = pa.struct([pa.field("protein_name", pa.string()), pa.field("peptide_count", pa.int32())])
PEPTIDE_COUNTS = pa.struct([pa.field("unique_sequences", pa.int32()), pa.field("total_sequences", pa.int32())])
PG_IDENTITY = ("pg_accessions", "grouped_runs", "label")
PG_SCHEMA = pa.schema([
    pa.field("pg_id", pa.int64(), nullable=False),
    pa.field("pg_accessions", pa.list_(pa.string()), nullable=False),
    pa.field("pg_group_id", pa.int64()),
    pa.field("anchor_protein", pa.string(), nullable=False),
    pa.field("grouped_runs", pa.list_(pa.string()), nullable=False),
    pa.field("label", pa.string()),
    pa.field("intensity", pa.float32()),
    pa.field("is_decoy", pa.bool_(), nullable=False),
    pa.field("contaminant", pa.bool_()),
    pa.field("peptides", pa.list_(PEPTIDE_PER_PROTEIN), nullable=False),
    pa.field("peptide_counts", PEPTIDE_COUNTS),
])

# run.yaml
SAMPLE_CHANNEL = pa.struct([
    pa.field("sample_accession", pa.string()),
    pa.field("label", pa.string()),
    pa.field("biological_replicate", pa.int32(), nullable=True),
    pa.field("technical_replicate", pa.int32(), nullable=True),
    pa.field("channel_role", pa.string(), nullable=True),
    pa.field("mixture", pa.string(), nullable=True),
])
RUN_SCHEMA = pa.schema([
    pa.field("run_accession", pa.string(), nullable=False),
    pa.field("run_file_name", pa.string(), nullable=False),
    pa.field("file_name", pa.string()),
    pa.field("samples", pa.list_(SAMPLE_CHANNEL), nullable=False),
    pa.field("fraction", pa.string()),
    pa.field("enzymes", pa.list_(pa.string())),
    pa.field("acquisition_method", pa.string()),
])
LFQ_LABEL = "LFQ"

# sample.yaml: the required fields, the optional ones an SDRF can fill, and
# the design's condition and replicate as extra columns (extra_columns: true)
NOT_AVAILABLE = "not available"
SDRF_SAMPLE_FIELDS = {
    "organism": "characteristics[organism]",
    "organism_part": "characteristics[organism part]",
    "disease": "characteristics[disease]",
    "cell_line": "characteristics[cell line]",
    "cell_type": "characteristics[cell type]",
    "sex": "characteristics[sex]",
    "age": "characteristics[age]",
    "developmental_stage": "characteristics[developmental stage]",
    "ancestry": "characteristics[ancestry category]",
    "individual": "characteristics[individual]",
    "sample_description": "description",
}
SAMPLE_SCHEMA = pa.schema(
    [pa.field("sample_accession", pa.string(), nullable=False),
     pa.field("organism", pa.string(), nullable=False),
     pa.field("organism_part", pa.string(), nullable=False)]
    + [pa.field(name, pa.string()) for name in list(SDRF_SAMPLE_FIELDS)[2:]]
    + [pa.field("condition", pa.string()),
       pa.field("biological_replicate", pa.int32()),
       pa.field("source_name", pa.string())]
)

# What each TOPP tool of the workflow does, as (step_category, step_name)
STEPS = {
    "DecoyDatabase": ("database", "decoy_database"),
    "ProSE": ("database_search", "sequence_search"),
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
    ("ln(hyperscore)", "prose ln(hyperscore)", True),
    ("hyperscore_zscore", "prose hyperscore z-score", True),
    ("delta_score", "prose delta score", True),
    ("COMET:xcorr", "comet xcorr", True),
    ("COMET:deltaCn", "comet deltacn", True),
    ("expect", "comet expectation value", False),
    ("q-value", "q-value", False),
]


def mass_error_ppm(observed: float, calculated: Optional[float], charge: int, isotope_error=None) -> Optional[float]:
    """Precursor m/z error in ppm to the isotope peak the search engine matched.

    ``isotope_error`` is the engine's offset (Sage: -1 when the observed m/z is
    one isotope above the monoisotope); without it the error is to the monoisotope.
    """
    if not calculated:
        return None
    offset = int(isotope_error) * C13_DELTA / charge if isotope_error and charge else 0.0
    return 1e6 * (observed + offset - calculated) / calculated


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
            "mass_error_ppm": mass_error_ppm(observed, calculated, charge, _float(_meta(hit, "isotope_error"))),
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
            "output_views": {"IDFilter": ["psm"], "ProteomicsLFQ": ["pg"]}.get(tool, []),
        })
    database = next(
        (p["value"] for r in rows if r["tool_name"] in ("OpenMS/ProSE", "OpenMS/CometAdapter") for p in r["parameters"] if p["key"] == "database"),
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
        "output_views": ["psm", "provenance", "run"] + (["pg"] if params.get("analysis-mode", "LFQ") == "LFQ" else []),
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


# ----------------------------------------------------------------------------- protein groups and runs

def _stem(file_name: str) -> str:
    name = Path(str(file_name)).name
    return name[:-5] if name.lower().endswith(".mzml") else Path(name).stem


def _stripped(peptide: str) -> str:
    """A peptide of the MSstats table without its modifications (``.(Oxidation)``, ``[+16]``)."""
    return re.sub(r"\(.*?\)|\[.*?\]|[^A-Za-z]", "", str(peptide)).upper()


def _int(value) -> Optional[int]:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def pg_rows(msstats, fractions: Optional[Mapping[str, dict]] = None) -> List[dict]:
    """One QPX pg row per protein group and sample of ProteomicsLFQ's MSstats table.

    ``msstats`` is the table (ProteinName, PeptideSequence, Reference,
    Intensity); a group's intensity in a sample is the sum of its rows, as on
    the results pages. ``fractions`` (``fractions.fraction_samples``) puts the
    fractions of a sample into one quantification unit (``grouped_runs``).
    """
    import pandas as pd

    from src.common.contaminants import is_contaminant
    from src.common.fractions import sample_of
    from src.common.protein_table import DECOY_PREFIXES

    df = msstats.dropna(subset=["ProteinName", "Reference"]).copy()
    if df.empty:
        return []
    df["_sample"] = df["Reference"].map(lambda ref: sample_of(ref, fractions))
    df["_run"] = df["Reference"].map(_stem)
    df["_seq"] = df["PeptideSequence"].map(_stripped) if "PeptideSequence" in df else ""
    df["Intensity"] = pd.to_numeric(df["Intensity"], errors="coerce")
    rows = []
    for (group, _sample), grp in df.groupby(["ProteinName", "_sample"], sort=True):
        accessions = [a for a in dict.fromkeys(str(group).split(";")) if a]
        if not accessions:
            continue
        runs = sorted(grp["_run"].unique().tolist())
        total = grp["PeptideSequence"].nunique() if "PeptideSequence" in grp else 0
        unique = grp["_seq"].nunique() if "PeptideSequence" in grp else 0
        intensity = grp["Intensity"].sum(min_count=1)
        rows.append({
            "pg_id": derive_id([sorted(accessions), runs, LFQ_LABEL]),
            "pg_accessions": accessions,
            "pg_group_id": derive_id(sorted(accessions)),
            "anchor_protein": accessions[0],
            "grouped_runs": runs,
            "label": LFQ_LABEL,
            "intensity": None if intensity != intensity else float(intensity),
            "is_decoy": all(a.startswith(DECOY_PREFIXES) for a in accessions),
            "contaminant": bool(is_contaminant(";".join(accessions))),
            "peptides": [{"protein_name": a, "peptide_count": int(unique)} for a in accessions],
            "peptide_counts": {"unique_sequences": int(unique), "total_sequences": int(total)},
        })
    return rows


ISOBARIC_PROTEIN_CSV = "openms_design_protein_openms.csv"
_ABUNDANCE = re.compile(r"^abundance_sample(\d+)\[")


def isobaric_label(tag: str, label: str) -> str:
    """QPX's label of a reporter channel: ``TMT127N``, ``iTRAQ114``."""
    return f"iTRAQ{tag}" if str(label).startswith("itraq") else f"TMT{tag}"


def isobaric_channels(columns, label: str) -> List[tuple]:
    """(column, channel number, QPX label) of each abundance column of
    ProteinQuantifier's protein table, in channel order."""
    from src.common.isobaric import ISOBARIC_CHANNELS

    tags = ISOBARIC_CHANNELS.get(label, [])
    out = []
    for col in columns:
        m = _ABUNDANCE.match(str(col))
        if m:
            n = int(m.group(1))
            tag = tags[n - 1] if n <= len(tags) else str(col).rsplit("[", 1)[1].rstrip("]")
            out.append((col, n, isobaric_label(tag, label)))
    return sorted(out, key=lambda c: c[1])


def isobaric_pg_rows(proteins, label: str, runs: List[str]) -> List[dict]:
    """One QPX pg row per protein group and reporter channel of an isobaric run.

    ``proteins`` is ProteinQuantifier's protein table (protein, n_proteins,
    protein_score, n_peptides, ``abundance_sampleN[...]``); the accessions of a
    group are separated by ``/``. Every channel's quantity covers all files of
    the labelled set, so ``grouped_runs`` names them all and ``label`` tells the
    channels apart.
    """
    import pandas as pd

    from src.common.contaminants import is_contaminant
    from src.common.protein_table import DECOY_PREFIXES

    if proteins is None or proteins.empty or "protein" not in proteins:
        return []
    channels = isobaric_channels(proteins.columns, label)
    runs = sorted(dict.fromkeys(runs))
    rows = []
    for _, protein in proteins.iterrows():
        accessions = [a for a in dict.fromkeys(str(protein["protein"]).split("/")) if a and a != "nan"]
        if not accessions:
            continue
        n_peptides = _int(protein.get("n_peptides")) or 0
        for col, _n, channel in channels:
            value = pd.to_numeric(protein[col], errors="coerce")
            rows.append({
                "pg_id": derive_id([sorted(accessions), runs, channel]),
                "pg_accessions": accessions,
                "pg_group_id": derive_id(sorted(accessions)),
                "anchor_protein": accessions[0],
                "grouped_runs": runs,
                "label": channel,
                "intensity": None if value != value else float(value),
                "is_decoy": all(a.startswith(DECOY_PREFIXES) for a in accessions),
                "contaminant": bool(is_contaminant(";".join(accessions))),
                "peptides": [{"protein_name": a, "peptide_count": n_peptides} for a in accessions],
                "peptide_counts": {"unique_sequences": n_peptides, "total_sequences": n_peptides},
            })
    return rows


def read_isobaric_proteins(results_dir):
    """ProteinQuantifier's protein table of an isobaric run (None without one)."""
    import pandas as pd

    path = Path(results_dir, "quant_results", ISOBARIC_PROTEIN_CSV)
    if not path.is_file():
        return None
    try:
        return pd.read_csv(path, sep="\t", comment="#")
    except (OSError, ValueError):
        return None


def _search_enzyme(params: Mapping, provenance: Optional[List[dict]] = None) -> Optional[str]:
    """The enzyme the run searched with: the search step of its provenance (every
    effective parameter, defaults included), else its stored settings."""
    for row in provenance or []:
        if row.get("step_name") != "sequence_search":
            continue
        values = {p["key"]: p["value"] for p in row.get("parameters") or []}
        if row.get("config"):
            values = {**json.loads(row["config"]), **values}
        for key in ("Search:enzyme", "enzyme"):
            if values.get(key):
                return str(values[key])
    for source in (params, params.get("_defaults") or {}):
        for key, value in source.items():
            if key.endswith(":1:Search:enzyme") or key.endswith(":1:enzyme"):
                return str(value)
    return None


def run_rows(params: Mapping, msstats=None, fractions: Optional[Mapping[str, dict]] = None,
             provenance: Optional[List[dict]] = None, channels: Optional[List[tuple]] = None) -> List[dict]:
    """One QPX run row per input file, with its sample, replicate and fraction.

    ``channels`` (``isobaric_channels``) gives an isobaric run's files one
    sample per reporter channel, as each file holds every channel.
    """
    from src.common.fractions import sample_of

    files = [Path(str(f)).name for f in params.get("mzML-files", []) or []]
    replicate = {}
    if msstats is not None and "Reference" in msstats:
        refs = msstats.drop_duplicates("Reference")
        reps = refs["BioReplicate"] if "BioReplicate" in refs else [None] * len(refs)
        for ref, rep in zip(refs["Reference"], reps):
            name = Path(str(ref)).name
            replicate[_stem(name)] = _int(rep)
            if not any(_stem(f) == _stem(name) for f in files):
                files.append(name)
    fraction_of = {}
    if fractions:
        # The design's fraction number, by the order of a sample's files
        by_sample = {}
        for f in sorted(fractions):
            by_sample.setdefault(fractions[f]["sample"], []).append(f)
        for members in by_sample.values():
            for i, f in enumerate(members, start=1):
                fraction_of[_stem(f)] = str(i)
    enzyme = _search_enzyme(params, provenance)
    rows = []
    for f in dict.fromkeys(files):
        stem = _stem(f)
        mzml = f if f.lower().endswith(".mzml") else f"{stem}.mzML"
        sample = sample_of(mzml, fractions)
        if channels:
            samples = [{
                "sample_accession": f"sample{n}",
                "label": channel,
                "biological_replicate": None,
                "technical_replicate": None,
                "channel_role": None,
                "mixture": None,
            } for _col, n, channel in channels]
        else:
            samples = [{
                "sample_accession": sample,
                "label": LFQ_LABEL if params.get("analysis-mode", "LFQ") == "LFQ" else None,
                "biological_replicate": replicate.get(stem),
                "technical_replicate": None,
                "channel_role": None,
                "mixture": None,
            }]
        rows.append({
            "run_accession": stem,
            "run_file_name": stem,
            "file_name": f,
            "samples": samples,
            "fraction": fraction_of.get(stem, "1"),
            "enzymes": [enzyme] if enzyme else None,
            "acquisition_method": "DDA",
        })
    return rows


_ORGANISM = re.compile(r"\bOS=(.+?)(?:\s+(?:OX|GN|PE|SV)=|$)")


def fasta_organism(path) -> Optional[str]:
    """The organism most entries of a UniProt-style FASTA name (``OS=``), else None."""
    counts = {}
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(">"):
                    m = _ORGANISM.search(line)
                    if m:
                        name = m.group(1).strip()
                        counts[name] = counts.get(name, 0) + 1
    except OSError:
        return None
    return max(counts, key=counts.get) if counts else None


def sdrf_samples(paths: Iterable) -> List[dict]:
    """Every row of the SDRF files as ``{column: first value}`` (lower-cased columns)."""
    rows = []
    for path in paths:
        try:
            with open(path, "r", encoding="utf-8-sig", newline="") as f:
                lines = [r for r in csv.reader(f, delimiter="\t") if any(c.strip() for c in r)]
        except OSError:
            continue
        if not lines:
            continue
        columns = [re.sub(r"\s+", " ", c.strip().lower()) for c in lines[0]]
        for line in lines[1:]:
            row = {}
            for column, value in zip(columns, line):
                if value.strip() and column not in row:
                    row[column] = value.strip()
            rows.append(row)
    return rows


def _sdrf_row(sdrf: List[dict], run_file_name: str, label: Optional[str]) -> Optional[dict]:
    """The SDRF row of a run's file (and, for a labelled run, of its channel)."""
    for row in sdrf:
        if _stem(Path(row.get("comment[data file]", "")).name) != run_file_name:
            continue
        if label and label != LFQ_LABEL and row.get("comment[label]", "").replace(" ", "").upper() != label.upper():
            continue
        return row
    return None


def sample_rows(run_table: List[dict], params: Mapping, fractions: Optional[Mapping[str, dict]] = None,
                sdrf: Optional[List[dict]] = None, organism: Optional[str] = None) -> List[dict]:
    """One QPX sample row per sample of the run view.

    The sample's condition is the group set for it (``mzML-group-<file>``, or
    ``TMT-group-<channel>`` of an isobaric run), the SDRF row of its file and
    channel fills the characteristics, and ``organism`` (the database's) stands
    in where no SDRF names one. QPX requires an organism and organism part, so
    those are "not available" when nothing names them.
    """
    sdrf = sdrf or []
    out = {}
    for run in run_table:
        mzml = run.get("file_name") or f"{run['run_file_name']}.mzML"
        if not str(mzml).lower().endswith(".mzml"):
            mzml = f"{run['run_file_name']}.mzML"
        for sample in run["samples"]:
            accession = sample["sample_accession"]
            label = sample.get("label")
            row = _sdrf_row(sdrf, run["run_file_name"], label)
            if accession in out:
                # A sample's later fractions fill what its first file's SDRF row left open
                if row and not out[accession]["source_name"]:
                    out[accession].update(_characteristics(row, organism))
                continue
            if label and label != LFQ_LABEL:
                condition = params.get(f"TMT-group-{accession}")
            elif fractions and mzml in fractions:
                condition = fractions[mzml].get("condition")
            else:
                condition = params.get(f"mzML-group-{mzml}")
            out[accession] = {
                "sample_accession": accession,
                **_characteristics(row, organism),
                "condition": str(condition).strip() or None if condition is not None else None,
                "biological_replicate": sample.get("biological_replicate"),
            }
    return list(out.values())


def _characteristics(row: Optional[dict], organism: Optional[str]) -> dict:
    row = row or {}
    values = {name: row.get(column) for name, column in SDRF_SAMPLE_FIELDS.items()}
    for name, value in values.items():
        if value and value.lower() in ("not applicable", "not available"):
            values[name] = None
    values["organism"] = values["organism"] or organism or NOT_AVAILABLE
    values["organism_part"] = values["organism_part"] or NOT_AVAILABLE
    values["source_name"] = row.get("source name")
    return values


def _sdrf_files(workflow_dir) -> List[Path]:
    """The SDRF files imported with the workflow's data (PRIDE import)."""
    if workflow_dir is None:
        return []
    folder = Path(workflow_dir, "pride-import", "sdrf")
    return sorted(folder.glob("*")) if folder.is_dir() else []


def _workflow_dir(results_dir) -> Optional[Path]:
    """The workflow folder of a run's results (``<workflow>/runs/<id>/results``)."""
    path = Path(results_dir).resolve()
    if path.parent.parent.name == "runs":
        return path.parent.parent.parent
    return path.parent


def read_msstats(results_dir):
    """ProteomicsLFQ's MSstats table of a run (None without one)."""
    import pandas as pd

    csvs = sorted(Path(results_dir, "quant_results").glob("*.csv"))
    if not csvs:
        return None
    try:
        return pd.read_csv(csvs[0])
    except (OSError, ValueError):
        return None


# ----------------------------------------------------------------------------- export

def qpx_paths(results_dir, name: str = "run") -> dict:
    folder = Path(results_dir, QPX_DIR)
    return {"psm": folder / f"{name}.psm.parquet", "provenance": folder / f"{name}.provenance.parquet"}


def view_paths(results_dir, name: str = "run") -> dict:
    """The protein group and run views, which a run written before them lacks."""
    folder = Path(results_dir, QPX_DIR)
    return {"pg": folder / f"{name}.pg.parquet", "run": folder / f"{name}.run.parquet",
            "sample": folder / f"{name}.sample.parquet"}


def export_run(results_dir, workflow_dir, params: Mapping, log_text: str, app_version: str,
               to_proforma=None, name: str = "run") -> dict:
    """Write a run's QPX psm and provenance views. Returns the paths and the PSM count."""
    if to_proforma is None:
        from src.common.peptide_format import to_proforma
    paths = qpx_paths(results_dir, name)
    idxmls = sorted(Path(results_dir, "psm_filter").glob("*.idXML"))
    mzml_names = [Path(str(f)).name for f in params.get("mzML-files", []) or []]
    n = write_psms(idxmls, mzml_names, paths["psm"], to_proforma)
    steps = provenance_rows(log_text, workflow_dir, params, app_version)
    write_table(steps, PROVENANCE_SCHEMA, paths["provenance"], footer("provenance_file"))
    return {**paths, "psms": n, **export_views(results_dir, params, name, steps)}


def export_views(results_dir, params: Mapping, name: str = "run", provenance: Optional[List[dict]] = None) -> dict:
    """Write a run's QPX run and sample views and its pg view.

    ``provenance`` (the run's provenance rows, read from its provenance view
    when not given) names the search's enzyme.
    """
    from src.common.fractions import fraction_samples

    paths = view_paths(results_dir, name)
    if provenance is None:
        stored = qpx_paths(results_dir, name)["provenance"]
        provenance = read_provenance(stored) if stored.is_file() else []
    from src.common.isobaric import label_type

    lfq = params.get("analysis-mode", "LFQ") == "LFQ"
    msstats = read_msstats(results_dir) if lfq else None
    proteins = None if lfq else read_isobaric_proteins(results_dir)
    label = label_type(dict(params))
    channels = isobaric_channels(proteins.columns, label) if proteins is not None else None
    fractions = fraction_samples(results_dir, params)
    run_table = run_rows(params, msstats, fractions, provenance, channels)
    out = {"run": write_table(run_table, RUN_SCHEMA, paths["run"],
                              footer("run_file", {"primary_key": "run_file_name"}))}
    database = params.get("fasta-file")
    organism = fasta_organism(database) if database and Path(str(database)).is_file() else None
    samples = sample_rows(run_table, params, fractions, sdrf_samples(_sdrf_files(_workflow_dir(results_dir))), organism)
    out["sample"] = write_table(samples, SAMPLE_SCHEMA, paths["sample"],
                                footer("sample_file", {"primary_key": "sample_accession"}))
    rows = None
    if msstats is not None:
        rows = pg_rows(msstats, fractions)
    elif proteins is not None:
        rows = isobaric_pg_rows(proteins, label, [r["run_file_name"] for r in run_table])
    if rows is not None:
        out["pg"] = write_table(rows, PG_SCHEMA, paths["pg"], footer("pg_file", {
            "primary_key": "pg_id",
            "identity_composite": ",".join(PG_IDENTITY),
        }))
    return {**out, "protein_groups": len(rows or [])}


def views_missing(results_dir, params: Mapping, name: str = "run") -> bool:
    """Whether a run lacks a QPX view it has the results for: one written
    before the run or sample view, or an isobaric run written before its pg view."""
    paths = view_paths(results_dir, name)
    if not paths["run"].is_file() or not paths["sample"].is_file():
        return True
    if params.get("analysis-mode", "LFQ") != "LFQ" and not paths["pg"].is_file():
        return Path(results_dir, "quant_results", ISOBARIC_PROTEIN_CSV).is_file()
    return False


def app_version() -> str:
    """The app's version from settings.json (empty when it cannot be read)."""
    try:
        with open("settings.json", "r", encoding="utf-8") as f:
            return str(json.load(f).get("version", ""))
    except (OSError, ValueError):
        return ""


def read_provenance(path) -> List[dict]:
    return pq.read_table(str(path)).to_pylist()
