"""Single run of the Viewer: one mzML file of an experiment, read into the library's columns."""
from pathlib import Path

import pytest

from src.viewer import run_files
from src.viewer.catalog import CATALOG, FILE_CATALOG, FILE_COMPONENTS, FILE_PRESETS


def test_file_ids_name_the_run_and_the_file():
    member = run_files.file_id("run-1", "sample_A")
    assert run_files.is_file(member) and not run_files.is_file("run-1")
    assert run_files.split(member) == ("run-1", "sample_A")
    assert run_files.stems({"mzML-files": ["/x/a.mzML", "/y/b.mzML"]}) == ["a", "b"]
    assert run_files.stems({"mzML-files": "/x/a.mzML"}) == ["a"]


def test_single_run_components_are_apart_from_the_comparison():
    names = {c["name"] for c in FILE_CATALOG}
    assert names == set(FILE_COMPONENTS)
    assert not names & {c["name"] for c in CATALOG}
    assert {n for p in FILE_PRESETS for row in p["layout"] for n in row} <= names


def test_read_tables_gives_scans_peaks_peak_map_and_chromatogram():
    pytest.importorskip("pyopenms")
    path = min(Path("example-data/mzML").glob("*.mzML"), key=lambda p: p.stat().st_size)
    tables = run_files.read_tables(path)
    scans, peaks = tables["scans"], tables["peaks"]
    assert scans.height > 0 and set(scans["MSLevel"].unique()) >= {1, 2}
    assert scans["#Peaks"].sum() == peaks.height
    assert tables["ms1"].height == scans.filter(scans["MSLevel"] == 1)["#Peaks"].sum()
    chromatogram = tables["chromatogram"]
    assert set(chromatogram["series"].unique()) == {"TIC", "Base peak"}
    assert chromatogram.height == 2 * scans.filter(scans["MSLevel"] == 1).height
