"""The files of a run that Export results copies or downloads."""
import io
import zipfile

from src.common import result_export as rx
from src.workflow.RunManager import RunManager


def _run(runs, name):
    run_id = runs.new_run(name)
    results = runs.results_dir(run_id)
    (results / "quant_results").mkdir()
    (results / "quant_results" / "msstats.csv").write_text("a,b\n")
    (results / "insight_cache").mkdir()
    (results / "insight_cache" / "x.parquet").write_text("cache")
    return run_id


def test_run_files_leave_out_caches(tmp_path):
    runs = RunManager(tmp_path)
    run_id = _run(runs, "A")
    rels = {rel.as_posix() for _, rel in rx.run_files(runs, run_id)}
    assert "results/quant_results/msstats.csv" in rels
    assert "run.json" in rels
    assert not any("insight_cache" in r for r in rels)


def test_all_runs_zip_has_a_folder_per_run(tmp_path):
    runs = RunManager(tmp_path)
    a = _run(runs, "Same: name")
    b = _run(runs, "Same: name")
    names = zipfile.ZipFile(io.BytesIO(rx.zip_runs(runs, [a, b]))).namelist()
    assert "Same - name/results/quant_results/msstats.csv" in names
    assert "Same - name (2)/results/quant_results/msstats.csv" in names
    assert not any("insight_cache" in n for n in names)


def test_safe_name():
    assert rx.safe_name('a/b:c*?') == "a_b_c__"
    assert rx.safe_name("  ") == "run"
