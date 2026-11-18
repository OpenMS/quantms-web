"""One PSM table over every file of a run.

The workflow caches a PSM table per input file (``table_<stem>``). A study of
many files is easier to browse in one table: ``all_files_table_id`` puts the
cached rows of every file together, with the file each PSM comes from, and
caches that as one more table. Each row keeps its file's identifiers, so the
spectrum, sequence and fragment views of that file can show the picked PSM
(``all_files_hit``).
"""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import polars as pl

from src.common.peptide_format import format_sequence

# The merged table's own row id, picked under this identifier
ROW_FIELD = "psm_row"
SELECTION = "psm_all"
FILE_FIELD = "file_name"
SOURCE_FIELD = "source"
DISPLAY_FIELD = "sequence_display"
CLEAN_FIELD = "not_contaminant"
PREFIX = "psms_all2_"  # change when the table's layout changes
SEARCH_FIELDS = (FILE_FIELD, "sequence", DISPLAY_FIELD, "protein_accession")


def file_label(idxml: Path) -> str:
    """The input a filter result belongs to: its name without the step's suffix."""
    stem = Path(idxml).stem
    return stem[: -len("_filter")] if stem.endswith("_filter") else stem


def _tables(cache_dir: Path, idxml_files: list) -> list[tuple[Path, Path]]:
    """(idXML, cached table folder) of every file whose table is cached."""
    found = []
    for idxml in idxml_files:
        base = Path(cache_dir) / f"table_{Path(idxml).stem}"
        if (base / "manifest.json").is_file() and (base / "preprocessed" / "data.parquet").is_file():
            found.append((Path(idxml), base))
    return found


def merged_rows(cache_dir: Path, idxml_files: list) -> pl.DataFrame:
    """Every file's cached PSM rows, with the file name, source and a row id of their own."""
    frames = []
    for idxml, base in _tables(cache_dir, idxml_files):
        data = pl.read_parquet(base / "preprocessed" / "data.parquet")
        frames.append(data.with_columns(
            pl.lit(file_label(idxml)).alias(FILE_FIELD), pl.lit(idxml.stem).alias(SOURCE_FIELD),
        ))
    if not frames:
        return pl.DataFrame()
    merged = pl.concat(frames, how="diagonal_relaxed")
    return merged.with_row_index(ROW_FIELD).with_columns(pl.col(ROW_FIELD).cast(pl.Int64))


def all_files_table_id(cache_dir: Path, idxml_files: list, mode: str = "OpenMS") -> str | None:
    """Cache id of the table of every file's PSMs, built on first use; None without any cached table.

    A new run of any file (a new cache generation) or another sequence form
    builds it again.
    """
    from openms_insight import Table

    tables = _tables(cache_dir, idxml_files)
    if not tables:
        return None
    generations = [
        (idxml.stem, str(json.loads((base / "manifest.json").read_text()).get("cache_generation", "")))
        for idxml, base in tables
    ]
    digest = hashlib.sha1(json.dumps([mode, generations]).encode()).hexdigest()[:12]
    table_id = f"{PREFIX}{digest}"
    if (Path(cache_dir) / table_id).is_dir():
        return table_id
    for stale in Path(cache_dir).glob(f"{PREFIX}*"):
        shutil.rmtree(stale, ignore_errors=True)

    data = merged_rows(cache_dir, idxml_files)
    config = json.loads((tables[0][1] / "manifest.json").read_text())["config"]
    columns = list(config["column_definitions"])
    if mode != "OpenMS":
        data = data.with_columns(
            pl.col("sequence").map_elements(lambda s: format_sequence(s, mode), return_dtype=pl.Utf8)
            .alias(DISPLAY_FIELD)
        )
        columns = [{**c, "field": DISPLAY_FIELD} if c.get("field") == "sequence" else c for c in columns]
    # Rounded for display, as the one-file table, so the File column fits
    # without folding other columns into the rows
    rounding = {"mz": 4, "rt": 1}
    columns = [
        {**c, "formatter": "money", "formatterParams": {"thousand": False, "precision": rounding[c["field"]]}}
        if c.get("field") in rounding and "formatter" not in c else c
        for c in columns
    ]
    columns.insert(0, {"field": FILE_FIELD, "title": "File", "headerTooltip": "The input file of this PSM"})
    data = data.with_columns(
        (~pl.col("contaminant")).alias(CLEAN_FIELD) if "contaminant" in data.columns
        else pl.lit(True).alias(CLEAN_FIELD)
    )
    Table(
        cache_id=table_id,
        data=data.lazy(),
        cache_path=str(cache_dir),
        interactivity={SELECTION: ROW_FIELD},
        column_definitions=columns,
        initial_sort=config.get("initial_sort"),
        index_field=ROW_FIELD,
        layout=config.get("layout", "fitDataFill"),
        page_size=config.get("page_size", 100),
        toggle_filters=[{"label": "Hide contaminants", "column": CLEAN_FIELD, "default": False}],
        search=[c for c in SEARCH_FIELDS if c in data.columns],
        search_placeholder="File, peptide or protein",
    )
    return table_id


def all_files_hit(cache_dir: Path, table_id: str, row) -> dict | None:
    """The source idXML stem and identifiers (id_idx, file_index, scan_id) of merged row ``row``."""
    if row is None:
        return None
    data = pl.scan_parquet(Path(cache_dir) / table_id / "preprocessed" / "data.parquet")
    hit = data.filter(pl.col(ROW_FIELD) == int(row)).select(
        SOURCE_FIELD, "id_idx", "file_index", "scan_id",
    ).collect()
    if hit.is_empty():
        return None
    return hit.row(0, named=True)
