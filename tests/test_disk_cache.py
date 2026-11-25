"""Viewer results kept on disk with their run (src/common/disk_cache.py)."""
import json
import os
import time

from src.common import disk_cache


def test_cache_root_is_the_run_folder(tmp_path):
    run = tmp_path / "wf" / "runs" / "r1"
    (run / "results" / "psm").mkdir(parents=True)
    (run / "run.json").write_text("{}")
    (tmp_path / "wf" / "params.json").write_text("{}")
    assert disk_cache.cache_root(run / "results" / "psm" / "a.idXML") == run
    other = tmp_path / "wf" / "input-files" / "a.mzML"
    other.parent.mkdir(parents=True)
    assert disk_cache.cache_root(other) == tmp_path / "wf"


def test_per_file_computes_once_and_again_when_the_file_changes(tmp_path):
    (tmp_path / "run.json").write_text("{}")
    source = tmp_path / "a.txt"
    source.write_text("one")
    calls = []

    @disk_cache.per_file
    def length(path, extra=0):
        calls.append(path)
        return len(open(path).read()) + extra

    assert length(source) == 3
    assert length(source) == 3
    assert len(calls) == 1
    assert length(source, 1) == 4
    assert len(calls) == 2
    assert list((tmp_path / disk_cache.FOLDER).glob("length-*.pkl"))

    time.sleep(0.01)
    source.write_text("three")
    os.utime(source)
    assert length(source) == 5
    assert len(calls) == 3


def test_unreadable_cache_is_computed_again(tmp_path):
    target = disk_cache.remember(tmp_path, "x", 1, lambda: {"a": 1})
    assert target == {"a": 1}
    for path in (tmp_path / disk_cache.FOLDER).glob("x-*.pkl"):
        path.write_bytes(b"not a pickle")
    assert disk_cache.remember(tmp_path, "x", 1, lambda: json.loads('{"b": 2}')) == {"b": 2}


def test_missing_file_is_not_cached(tmp_path):
    @disk_cache.per_file
    def name(path):
        return str(path)

    assert name(tmp_path / "missing") == str(tmp_path / "missing")
    assert not (tmp_path / disk_cache.FOLDER).exists()
