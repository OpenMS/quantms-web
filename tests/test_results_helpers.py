from pathlib import Path

from src.common.results_helpers import (
    extract_filename_from_idxml,
    extract_scan_from_ref,
    extract_scan_number,
    get_workflow_dir,
)


def test_get_workflow_dir():
    assert get_workflow_dir("/tmp/ws") == Path("/tmp/ws", "topp-workflow")


def test_extract_scan_from_ref():
    assert extract_scan_from_ref("controllerType=0 controllerNumber=1 scan=1234") == 1234
    assert extract_scan_from_ref("no scan token here") == 0


def test_extract_scan_number():
    assert extract_scan_number("file=abc scan=42") == 42


def test_extract_filename_from_idxml_strips_suffixes():
    assert extract_filename_from_idxml(Path("02COVID_filter.idXML")) == "02COVID.mzML"
    assert extract_filename_from_idxml(Path("sample_comet.idXML")) == "sample.mzML"
    assert extract_filename_from_idxml(Path("run_per.idXML")) == "run.mzML"
