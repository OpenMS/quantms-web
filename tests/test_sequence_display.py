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


def test_filtered_table_adds_q_values_and_searches(tmp_path):
    import pyopenms as poms
    from src.common.sequence_display import CLEAN_FIELD, Q_FIELD, Q_IDENTIFIER, filtered_table_id

    _base_table(tmp_path)
    # The run's idXML, hits in id_idx order with their q-values
    peptides = []
    for seq, q in (("PEPM(Oxidation)K", 0.004), ("PEPTIDEK", 0.0002)):
        pep = poms.PeptideIdentification()
        hit = poms.PeptideHit()
        hit.setSequence(poms.AASequence.fromString(seq))
        hit.setMetaValue("MS:1001491", q)
        pep.setHits([hit])
        peptides.append(pep)
    if hasattr(poms, "PeptideIdentificationList"):
        container = poms.PeptideIdentificationList()
        for pep in peptides:
            container.push_back(pep)
        peptides = container
    idxml = tmp_path / "run.idXML"
    poms.IdXMLFile().store(str(idxml), [poms.ProteinIdentification()], peptides)

    variant = filtered_table_id(tmp_path, "table_run", "ProForma", idxml)
    data = pl.read_parquet(tmp_path / variant / "preprocessed" / "data.parquet").sort("id_idx")
    assert data[Q_FIELD].to_list() == pytest.approx([0.004, 0.0002])
    assert data[CLEAN_FIELD].to_list() == [True, True]
    manifest = json.loads((tmp_path / variant / "manifest.json").read_text())
    assert manifest["config"]["range_filters"] == {Q_IDENTIFIER: [Q_FIELD, "q_ceiling"]}
    assert manifest["config"]["toggle_filters"][0]["column"] == CLEAN_FIELD
    fields = [c["field"] for c in manifest["config"]["column_definitions"]]
    assert fields.index(Q_FIELD) == fields.index("score") + 1

    # A search keeps the matching rows (any modification form) and replaces the variant
    searched = filtered_table_id(tmp_path, "table_run", "ProForma", idxml, "m[oxidation]")
    assert searched != variant and not (tmp_path / variant).exists()
    data = pl.read_parquet(tmp_path / searched / "preprocessed" / "data.parquet")
    assert data["id_idx"].to_list() == [0]
