"""PSM table with a switch for how modified peptides are written.

The workflow caches each run's PSM table with OpenMS sequence strings. When the
user picks another form (see src/common/peptide_format.py), the cached table's
data is reused for a second table cache with that form in the Sequence column.
It is built once per run and form, keeps the same row ids and links, so the
spectrum and sequence views follow it as they follow the original.
"""
import json
import shutil
from pathlib import Path

import polars as pl
import streamlit as st
from openms_insight import Table

from src.common.peptide_format import DISPLAY_MODES, format_sequence

# The choice lives under MODE_KEY and is shared by every PSM table. Each table
# has its own radio (several can be on one page), seeded from MODE_KEY on each
# render, so changing one changes them all.
MODE_KEY = "sequence-display-mode"
DISPLAY_FIELD = "sequence_display"


def sequence_mode_selector(name: str) -> str:
    """Radio for the modification display of the PSM table ``name``."""
    from src.viewer.scope import scoped

    widget_key = scoped(f"_{MODE_KEY}-{name}")
    st.session_state[widget_key] = st.session_state.get(MODE_KEY, DISPLAY_MODES[0])

    def _changed() -> None:
        st.session_state[MODE_KEY] = st.session_state[widget_key]

    mode = st.radio(
        "Modifications",
        DISPLAY_MODES,
        key=widget_key,
        on_change=_changed,
        horizontal=True,
        help=(
            "Labels: modified residues carry a short label, e.g. M(ox), S(ph), C(cam). "
            "ProForma: the HUPO-PSI standard, e.g. M[Oxidation]; PSM downloads always include it. "
            "OpenMS: the string as OpenMS writes it."
        ),
    )
    st.session_state[MODE_KEY] = mode
    return mode


def display_table_id(cache_dir: Path, table_id: str, mode: str) -> str:
    """Cache id of ``table_id`` with its sequences in ``mode``, built on first use."""
    if mode == "OpenMS":
        return table_id
    base = Path(cache_dir) / table_id
    manifest = json.loads((base / "manifest.json").read_text())
    # A new run regenerates the base cache, which gets a new generation id
    generation = str(manifest.get("cache_generation", ""))[:8]
    variant_prefix = f"{table_id}__{mode.lower()}_"
    variant_id = f"{variant_prefix}{generation}"
    if (Path(cache_dir) / variant_id).is_dir():
        return variant_id

    for stale in Path(cache_dir).glob(f"{variant_prefix}*"):
        shutil.rmtree(stale, ignore_errors=True)
    data = pl.read_parquet(base / "preprocessed" / "data.parquet")
    data = data.with_columns(
        pl.col("sequence")
        .map_elements(lambda s: format_sequence(s, mode), return_dtype=pl.Utf8)
        .alias(DISPLAY_FIELD)
    )
    config = manifest["config"]
    column_definitions = [
        {**col, "field": DISPLAY_FIELD} if col.get("field") == "sequence" else col
        for col in config["column_definitions"]
    ]
    Table(
        cache_id=variant_id,
        data=data.lazy(),
        cache_path=str(cache_dir),
        interactivity=manifest.get("interactivity") or None,
        column_definitions=column_definitions,
        initial_sort=config.get("initial_sort"),
        index_field=config.get("index_field", "id"),
        go_to_fields=config.get("go_to_fields"),
        layout=config.get("layout", "fitDataFill"),
        page_size=config.get("page_size", 100),
    )
    return variant_id


def psm_table(
    cache_dir: Path, table_id: str, mode: str, state_manager, height: int = 533, title: str | None = "PSMs",
) -> None:
    """Render the cached PSM table ``table_id`` with its sequences in ``mode``."""
    table = Table(cache_id=display_table_id(cache_dir, table_id, mode), cache_path=str(cache_dir), collapse_columns=False)
    if title:
        table._title = title
    table(state_manager=state_manager, height=height)


# The PSM page's table menu: a q-value column with thresholds, a contaminant
# switch, and a search over the sequence (in every written form) and proteins
Q_FIELD = "q_value"
CLEAN_FIELD = "not_contaminant"
SEARCH_FIELDS = ("sequence", DISPLAY_FIELD, "proforma", "protein_accession")
SEARCH_PLACEHOLDER = "Peptide or protein"


def filtered_table_id(
    cache_dir: Path, table_id: str, mode: str, idxml_path: Path, q_suggestions: list | None = None,
) -> str:
    """Cache id of the PSM table with a q-value column, a contaminant switch, a
    search box and q-value thresholds in its filter dialog; built on first use,
    from the table in ``mode``.

    The table filters itself (in Python, over every row), so neither a search
    nor a threshold needs a new cache.
    """
    import hashlib

    from src.common.result_filters import q_value_per_hit
    from src.common.results_helpers import load_idxml

    source_id = display_table_id(cache_dir, table_id, mode)
    source = Path(cache_dir) / source_id
    manifest = json.loads((source / "manifest.json").read_text())
    generation = str(manifest.get("cache_generation", ""))[:8]
    q_suggestions = list(q_suggestions or [])
    digest = hashlib.sha1(f"{source_id}|{json.dumps(q_suggestions, sort_keys=True)}".encode()).hexdigest()[:8]
    prefix = f"{table_id}__f3"  # bump when the variant's layout changes
    variant_id = f"{prefix}{digest}_{generation}"
    if (Path(cache_dir) / variant_id).is_dir():
        return variant_id

    for stale in Path(cache_dir).glob(f"{table_id}__f*"):  # every older filtered variant
        shutil.rmtree(stale, ignore_errors=True)
    data = pl.read_parquet(source / "preprocessed" / "data.parquet")
    q_values = q_value_per_hit(load_idxml(str(idxml_path))[1])
    ids = data["id_idx"].to_numpy()
    in_range = (ids >= 0) & (ids < len(q_values))
    q = [float(q_values[i]) if ok else None for i, ok in zip(ids, in_range)]
    data = data.with_columns(
        pl.Series(Q_FIELD, q, dtype=pl.Float64).fill_nan(None),
        (~pl.col("contaminant")).alias(CLEAN_FIELD) if "contaminant" in data.columns else pl.lit(True).alias(CLEAN_FIELD),
    )

    config = manifest["config"]
    # Rounded for display (the cache holds float32 noise like 494.2902526855469),
    # which keeps the table narrow enough not to fold columns into the rows
    rounding = {"mz": 4, "rt": 1}
    column_definitions = [
        {**c, "formatter": "money", "formatterParams": {"thousand": False, "precision": rounding[c["field"]]}}
        if c.get("field") in rounding and "formatter" not in c else c
        for c in config["column_definitions"]
    ]
    q_column = {"field": Q_FIELD, "title": "q-value", "sorter": "number",
                "formatter": "money", "formatterParams": {"thousand": False, "precision": 5},
                "headerTooltip": "Percolator's q-value: the lowest FDR at which this match is accepted"}
    after = next((i for i, c in enumerate(column_definitions) if c.get("field") == "score"), len(column_definitions) - 1)
    column_definitions.insert(after + 1, q_column)
    Table(
        cache_id=variant_id,
        data=data.lazy(),
        cache_path=str(cache_dir),
        interactivity=manifest.get("interactivity") or None,
        column_definitions=column_definitions,
        initial_sort=config.get("initial_sort"),
        index_field=config.get("index_field", "id"),
        go_to_fields=config.get("go_to_fields"),
        layout=config.get("layout", "fitDataFill"),
        page_size=config.get("page_size", 100),
        toggle_filters=[{"label": "Hide contaminants", "column": CLEAN_FIELD, "default": False}],
        search=[c for c in SEARCH_FIELDS if c in data.columns],
        search_placeholder=SEARCH_PLACEHOLDER,
        numeric_filters={Q_FIELD: {"op": "<=", "suggestions": q_suggestions}},
    )
    return variant_id


def filtered_psm_table(
    cache_dir: Path, table_id: str, mode: str, state_manager, idxml_path: Path,
    q_suggestions: list | None = None, height: int = 533, title: str | None = "PSMs",
) -> None:
    """Render the PSM table with its search box, q-value thresholds and contaminant switch."""
    variant_id = filtered_table_id(cache_dir, table_id, mode, idxml_path, q_suggestions)
    table = Table(cache_id=variant_id, cache_path=str(cache_dir), collapse_columns=False)
    if title:
        table._title = title
    table(state_manager=state_manager, height=height)
