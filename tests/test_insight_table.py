"""Tests for the OpenMS-Insight table helpers in src/common/insight_tables.py."""
import numpy as np
import pandas as pd
import polars as pl

from src.common.insight_tables import ROW_FIELD, _fingerprint_dataframe, _to_polars


def _abundance_df():
    return pd.DataFrame(
        {
            "ProteinName": ["P1", "P2", "P3"],
            "m/z": [1.5, np.nan, 3.5],
            "sample1[126.1]": [1, 2, 3],
            "PeptideSequence": ["AAK", "BBK", "CCK"],
        }
    )


def test_to_polars_adds_unique_hidden_row_id():
    pl_df, defs = _to_polars(_abundance_df())
    assert pl_df[ROW_FIELD].n_unique() == 3
    assert ROW_FIELD not in [d["field"] for d in defs]


def test_to_polars_keeps_titles_and_sanitizes_dotted_fields():
    pl_df, defs = _to_polars(_abundance_df())
    by_title = {d["title"]: d for d in defs}
    # Tabulator reads "." in a field name as a nested path, so it is replaced;
    # the visible title keeps the original column name.
    assert by_title["sample1[126.1]"]["field"] == "sample1[126_1]"
    assert by_title["m/z"]["field"] == "m/z"
    assert set(d["field"] for d in defs) <= set(pl_df.columns)


def test_to_polars_makes_duplicate_column_names_unique():
    df = pd.DataFrame([[1, 2]], columns=["a.b", "a_b"])
    _, defs = _to_polars(df)
    assert len({d["field"] for d in defs}) == 2


def test_to_polars_numeric_columns_sort_as_numbers():
    _, defs = _to_polars(_abundance_df())
    by_title = {d["title"]: d for d in defs}
    assert by_title["m/z"]["sorter"] == "number"
    assert "sorter" not in by_title["ProteinName"]


def test_to_polars_turns_named_index_into_column():
    df = _abundance_df().set_index("ProteinName")
    pl_df, defs = _to_polars(df)
    assert "ProteinName" in pl_df.columns
    assert defs[0]["title"] == "ProteinName"


def test_to_polars_drops_unnamed_index_after_filtering():
    pl_df, _ = _to_polars(_abundance_df().iloc[[0, 2]])
    assert "index" not in pl_df.columns


def test_to_polars_nan_becomes_null_and_mixed_objects_are_stringified():
    df = pd.DataFrame({"x": [1.0, np.nan], "mixed": ["a", 1]})
    pl_df, _ = _to_polars(df)
    assert pl_df["x"].null_count() == 1
    assert pl_df["mixed"].dtype == pl.String


def test_fingerprint_is_stable_and_detects_changes():
    df = _abundance_df()
    assert _fingerprint_dataframe(df) == _fingerprint_dataframe(df.copy())

    changed_value = df.copy()
    changed_value.loc[0, "sample1[126.1]"] = 99
    assert _fingerprint_dataframe(changed_value) != _fingerprint_dataframe(df)

    assert _fingerprint_dataframe(df.iloc[:2]) != _fingerprint_dataframe(df)
    assert _fingerprint_dataframe(df.rename(columns={"m/z": "mz"})) != _fingerprint_dataframe(df)


def test_fingerprint_handles_duplicate_columns_and_list_cells():
    _fingerprint_dataframe(pd.DataFrame([[1, 2]], columns=["a", "a"]))
    _fingerprint_dataframe(pd.DataFrame({"x": [[1, 2], [3]]}))
    _fingerprint_dataframe(pd.DataFrame({"x": []}))


def test_table_without_source_files_is_cached_by_content(tmp_path, monkeypatch):
    from src.common import insight_tables

    rendered = []
    monkeypatch.setattr(insight_tables.Table, "__call__", lambda self, **kw: rendered.append(self._cache_id))
    df = _abundance_df()
    insight_tables.show_insight_table(df, "t", tmp_path)
    insight_tables.show_insight_table(df.copy(), "t", tmp_path)
    assert rendered[0] == rendered[1]

    changed = df.copy()
    changed.loc[0, "sample1[126.1]"] = 99
    insight_tables.show_insight_table(changed, "t", tmp_path)
    assert rendered[2] != rendered[0]
    # The cache of the previous version of the table is removed
    assert [p.name for p in tmp_path.iterdir()] == [rendered[2]]
