"""The database search runs one ProSE process per share of the threads."""
from types import SimpleNamespace

# Imported inside the tests, as in test_contaminants.py: importing the
# workflow at collection time would give later tests that mock streamlit the
# real src.workflow modules.


def search_batches(n_files, max_threads):
    from src.WorkflowTest import search_batches as batches

    return batches(n_files, max_threads)


def test_one_process_below_two_shares_of_threads():
    assert search_batches(6, 8) == [[0, 1, 2, 3, 4, 5]]
    assert search_batches(200, 31) == [list(range(200))]


def test_files_dealt_out_evenly_and_in_order():
    batches = search_batches(200, 64)
    assert len(batches) == 4
    assert [len(b) for b in batches] == [50, 50, 50, 50]
    assert sum(batches, []) == list(range(200))
    assert [len(b) for b in search_batches(10, 48)] == [4, 3, 3]


def test_never_more_processes_than_files():
    assert search_batches(2, 128) == [[0], [1]]
    assert search_batches(0, 64) == []
    assert search_batches(3, 1) == [[0, 1, 2]]


def test_run_search_passes_each_batch_as_one_command(tmp_path):
    from src.WorkflowTest import WorkflowTest

    calls = []
    wf = object.__new__(WorkflowTest)
    wf.logger = SimpleNamespace(log=lambda *a, **k: None)
    wf.executor = SimpleNamespace(
        _get_max_threads=lambda: 32,
        run_topp=lambda tool, io, params, tool_instance_name=None: calls.append((tool, io, params)) or True,
    )
    mzml = [f"r{i}.mzML" for i in range(5)]
    out = [f"r{i}.idXML" for i in range(5)]

    assert wf.run_search("ProSE", mzml, out, tmp_path / "db.fasta")

    tool, io, params = calls[0]
    assert tool == "ProSE"
    assert io["in"] == [mzml[:3], mzml[3:]]
    assert io["out_idxml"] == [out[:3], out[3:]]
    assert params["database"] == str(tmp_path / "db.fasta")
