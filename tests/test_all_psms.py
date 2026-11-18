"""The PSM table over every file of a run (src/common/all_psms.py)."""
import pytest

from tests.test_identified_spectra import _write_idxml, _write_mzml


def _run(tmp_path, files):
    from src.common.identified_spectra import write_identified_spectra
    from src.common.psm_views import rebuild_psm_views

    inputs = tmp_path / "in"
    inputs.mkdir()
    results = tmp_path / "results"
    (results / "psm_filter").mkdir(parents=True)
    mzmls, idxmls = [], []
    for name, scans in files.items():
        mzml = inputs / f"{name}.mzML"
        _write_mzml(mzml)
        idxml = results / "psm_filter" / f"{name}_filter.idXML"
        _write_idxml(idxml, scans=scans)
        mzmls.append(mzml)
        idxmls.append(idxml)
    write_identified_spectra(mzmls, idxmls, results)
    for idxml in idxmls:
        assert rebuild_psm_views(results, idxml, mzmls, 0.02, False)
    return results / "insight_cache", idxmls


def test_one_table_holds_every_file(tmp_path):
    pytest.importorskip("openms_insight")
    import polars as pl

    from src.common.all_psms import all_files_hit, all_files_table_id, FILE_FIELD, ROW_FIELD

    cache, idxmls = _run(tmp_path, {"a": [3, 5], "b": [2, 3, 4]})
    table_id = all_files_table_id(cache, idxmls)
    data = pl.read_parquet(cache / table_id / "preprocessed" / "data.parquet")
    assert data.height == 5
    assert sorted(data[ROW_FIELD].to_list()) == list(range(5))
    assert data.group_by(FILE_FIELD).len().sort(FILE_FIELD)["len"].to_list() == [2, 3]
    # A picked row leads back to its own file's identifiers
    row = data.filter(pl.col(FILE_FIELD) == "b").sort("scan_id")[ROW_FIELD][0]
    hit = all_files_hit(cache, table_id, row)
    assert hit["source"] == "b_filter"
    assert hit["file_index"] == 1 and hit["scan_id"] == 2
    assert all_files_hit(cache, table_id, None) is None
    # Built once; another sequence form is a table of its own and replaces it
    assert all_files_table_id(cache, idxmls) == table_id
    labels = all_files_table_id(cache, idxmls, mode="Labels")
    assert labels != table_id
    assert not (cache / table_id).exists()
    assert "sequence_display" in pl.read_parquet(cache / labels / "preprocessed" / "data.parquet").columns


def test_files_without_a_table_are_left_out(tmp_path):
    pytest.importorskip("openms_insight")
    from src.common.all_psms import all_files_table_id, file_label

    cache, idxmls = _run(tmp_path, {"a": [3]})
    assert all_files_table_id(cache, [idxmls[0].with_name("c_filter.idXML")]) is None
    assert all_files_table_id(cache, [*idxmls, idxmls[0].with_name("c_filter.idXML")]) is not None
    assert file_label(idxmls[0]) == "a"
