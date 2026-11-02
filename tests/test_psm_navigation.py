"""Previous/Next through a PSM table in its current sort and filters."""
import polars as pl
import pytest

pytest.importorskip("openms_insight")

# src.common.psm_navigation imports streamlit, so it is imported inside the tests,
# after test_parameter_defaults has imported its modules with a mock.


def _table(tmp_path):
    from openms_insight import Table

    data = pl.DataFrame({
        "id_idx": [0, 1, 2, 3],
        "file_index": [0, 0, 0, 0],
        "scan_id": [10, 11, 12, 13],
        "peptide": ["A", "B", "A", "C"],
        "mz": [500.0, 400.0, 600.0, 450.0],
    })
    Table(
        cache_id="psms", data=data.lazy(), cache_path=str(tmp_path),
        filters={"peptide": "peptide"},
        interactivity={"identification": "id_idx", "file": "file_index", "spectrum": "scan_id"},
        index_field="id_idx", initial_sort=[{"column": "id_idx", "dir": "asc"}],
    )
    return Table(cache_id="psms", cache_path=str(tmp_path))


INTERACTIVITY = {"identification": "id_idx", "file": "file_index", "spectrum": "scan_id"}


def test_steps_follow_the_initial_order(tmp_path):
    from src.common.psm_navigation import neighbour, row_position

    table = _table(tmp_path)
    state = {"identification": 1}
    assert row_position(table, state, "id_idx", 1)[0] == 1
    assert neighbour(table, state, INTERACTIVITY, "identification", 1) == {
        "identification": 2, "file": 0, "spectrum": 12,
    }
    assert neighbour(table, state, INTERACTIVITY, "identification", -1)["identification"] == 0
    # At the ends there is no neighbour
    assert neighbour(table, {"identification": 3}, INTERACTIVITY, "identification", 1) is None
    assert neighbour(table, {"identification": 0}, INTERACTIVITY, "identification", -1) is None


def test_steps_follow_the_viewers_sort_and_filters(tmp_path):
    from src.common.psm_navigation import neighbour, row_position

    table = _table(tmp_path)
    # The viewer sorted by m/z descending: 600 (2), 500 (0), 450 (3), 400 (1)
    state = {"identification": 0, "psms_page": {"sort_column": "mz", "sort_dir": "desc"}}
    assert row_position(table, state, "id_idx", 0.0)[0] == 1  # the browser sends floats
    assert neighbour(table, state, INTERACTIVITY, "identification", 1)["identification"] == 3
    assert neighbour(table, state, INTERACTIVITY, "identification", -1)["identification"] == 2
    # A linked filter (the selected peptide) limits the rows stepped through
    state = {"identification": 0, "peptide": "A"}
    index, values = row_position(table, state, "id_idx", 0)
    assert (index, values) == (0, [0, 2])
    assert neighbour(table, state, INTERACTIVITY, "identification", 1)["identification"] == 2
    # A row the filters hide has no position
    assert row_position(table, {"peptide": "A"}, "id_idx", 1)[0] is None
