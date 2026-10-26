"""Tests for the PSM table variants with another modification display."""
import json

import polars as pl
import pytest

pytest.importorskip("openms_insight")

# src.common.sequence_display imports streamlit, so it is imported inside the
# tests, after test_parameter_defaults has imported its modules with a mock.


def _base_table(cache_dir):
    from openms_insight import Table

    data = pl.DataFrame({
        "id_idx": [0, 1],
        "scan_id": [1, 2],
        "file_index": [0, 0],
        "sequence": ["PEPM(Oxidation)K", "PEPTIDEK"],
        "score": [0.01, 0.02],
    })
    Table(
        cache_id="table_run",
        data=data.lazy(),
        cache_path=str(cache_dir),
        interactivity={"identification": "id_idx"},
        column_definitions=[
            {"field": "sequence", "title": "Sequence"},
            {"field": "score", "title": "Score", "sorter": "number"},
        ],
        initial_sort=[{"column": "score", "dir": "asc"}],
        index_field="id_idx",
    )


def test_variant_table_shows_the_chosen_form(tmp_path):
    from src.common.sequence_display import DISPLAY_FIELD, display_table_id

    _base_table(tmp_path)
    assert display_table_id(tmp_path, "table_run", "OpenMS") == "table_run"

    variant = display_table_id(tmp_path, "table_run", "ProForma")
    assert variant.startswith("table_run__proforma_")
    data = pl.read_parquet(tmp_path / variant / "preprocessed" / "data.parquet")
    assert data.sort("id_idx")[DISPLAY_FIELD].to_list() == ["PEPM[Oxidation]K", "PEPTIDEK"]
    manifest = json.loads((tmp_path / variant / "manifest.json").read_text())
    # Same links and row ids as the original, with the Sequence column swapped
    assert manifest["interactivity"] == {"identification": "id_idx"}
    assert manifest["config"]["index_field"] == "id_idx"
    assert manifest["config"]["column_definitions"][0] == {"field": DISPLAY_FIELD, "title": "Sequence"}

    # Built once: asking again reuses it
    assert display_table_id(tmp_path, "table_run", "ProForma") == variant


def test_variant_follows_a_rebuilt_base(tmp_path):
    from openms_insight import Table
    from src.common.sequence_display import display_table_id

    _base_table(tmp_path)
    first = display_table_id(tmp_path, "table_run", "Labels")
    # A new run regenerates the base table
    Table(
        cache_id="table_run",
        data=pl.DataFrame({"id_idx": [0], "sequence": ["ELVISK"], "score": [0.5]}).lazy(),
        cache_path=str(tmp_path),
        column_definitions=[{"field": "sequence", "title": "Sequence"}],
        index_field="id_idx",
        regenerate_cache=True,
    )
    second = display_table_id(tmp_path, "table_run", "Labels")
    assert second != first
    assert not (tmp_path / first).exists()
    data = pl.read_parquet(tmp_path / second / "preprocessed" / "data.parquet")
    assert data["sequence_display"].to_list() == ["ELVISK"]
