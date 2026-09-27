"""Streamed result tables built on openms-insight.

``st.dataframe`` ships the whole table to the browser, which stalls on
full-scale runs. An openms-insight ``Table`` caches the data on disk once and
sends one page at a time; its download button exports every row and column.
"""
import hashlib
import shutil
from pathlib import Path

import pandas as pd
import polars as pl
from openms_insight import Table

ROW_FIELD = "row"


def _cache_key(source_files, extra: str) -> str:
    """Hash the source files' identity and mtimes, so new results get a new cache."""
    parts = [extra]
    for path in source_files:
        path = Path(path)
        mtime = path.stat().st_mtime_ns if path.exists() else 0
        parts.append(f"{path}:{mtime}")
    return hashlib.sha1("|".join(parts).encode()).hexdigest()[:10]


def _to_polars(df: pd.DataFrame) -> tuple[pl.DataFrame, list[dict]]:
    """Convert to polars with Tabulator-safe field names and column definitions.

    Tabulator reads a "." in a field name as a nested lookup, so such columns get
    a safe field name and keep their original name as the column title.
    """
    df = df.reset_index(drop=True)
    renames, column_definitions, used = {}, [], {ROW_FIELD}
    for column in df.columns:
        field = str(column).replace(".", "_")
        while field in used:
            field = f"{field}_"
        used.add(field)
        renames[column] = field
        col_def = {"field": field, "title": str(column), "headerTooltip": True}
        if pd.api.types.is_numeric_dtype(df[column]):
            col_def["sorter"] = "number"
            col_def["hozAlign"] = "right"
        column_definitions.append(col_def)
    df = df.rename(columns=renames)
    for column in df.columns:
        if df[column].dtype == object:
            df[column] = df[column].map(lambda v: None if v is None or v != v else str(v))
    pl_df = pl.from_pandas(df).with_row_index(ROW_FIELD)
    return pl_df, column_definitions


def show_insight_table(
    df: pd.DataFrame,
    name: str,
    cache_dir,
    source_files=(),
    extra_key: str = "",
    title: str | None = None,
    height: int = 500,
    page_size: int = 100,
    **table_kwargs,
) -> None:
    """Render ``df`` as a server-side paginated openms-insight table.

    The cache is keyed on ``source_files`` (path and mtime) and ``extra_key``;
    when either changes (a rerun, new group assignments) a new cache replaces the
    old one, so the page never shows a previous run's rows.

    Args:
        df: Table to show.
        name: Stable name for this table; must be unique per page.
        cache_dir: Directory holding openms-insight caches.
        source_files: Files the table is derived from.
        extra_key: Anything else the table depends on (e.g. parameters).
        title: Title, also used as the CSV download's file name.
        height: Component height in pixels.
        page_size: Rows per page.
        **table_kwargs: Passed on to ``Table`` (e.g. ``initial_sort``).
    """
    cache_dir = Path(cache_dir)
    cache_id = f"{name}_{_cache_key(source_files, extra_key)}"

    if not (cache_dir / cache_id).is_dir():
        # Drop this table's caches from earlier results before building the new one
        for stale in cache_dir.glob(f"{name}_" + "?" * 10):
            shutil.rmtree(stale, ignore_errors=True)
        cache_dir.mkdir(parents=True, exist_ok=True)
        pl_df, column_definitions = _to_polars(df)
        Table(
            cache_id=cache_id,
            data=pl_df.lazy(),
            cache_path=str(cache_dir),
            column_definitions=column_definitions,
            index_field=ROW_FIELD,
            title=title or name,
            page_size=page_size,
            **table_kwargs,
        )

    table = Table(cache_id=cache_id, cache_path=str(cache_dir))
    table(key=f"insight_table_{name}", height=height)
