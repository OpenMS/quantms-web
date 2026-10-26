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

# The choice lives under MODE_KEY; the radio's own key is dropped by Streamlit
# whenever a page without it runs, so it is seeded from MODE_KEY on each render
MODE_KEY = "sequence-display-mode"
WIDGET_KEY = f"_{MODE_KEY}"
DISPLAY_FIELD = "sequence_display"


def sequence_mode_selector() -> str:
    """Radio for the modification display; the choice is shared by every PSM page."""
    if WIDGET_KEY not in st.session_state:
        st.session_state[WIDGET_KEY] = st.session_state.get(MODE_KEY, DISPLAY_MODES[0])
    mode = st.radio(
        "Modifications",
        DISPLAY_MODES,
        key=WIDGET_KEY,
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


def psm_table(cache_dir: Path, table_id: str, mode: str, state_manager, height: int = 533) -> None:
    """Render the cached PSM table ``table_id`` with its sequences in ``mode``."""
    table = Table(cache_id=display_table_id(cache_dir, table_id, mode), cache_path=str(cache_dir))
    table(state_manager=state_manager, height=height)
