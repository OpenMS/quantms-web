"""Keep only final results: the search and rescoring hits go, the final results stay."""
from pathlib import Path

from src.common.final_results import drop_intermediates, keeps_final_only, size_label


def _write(path: Path, size: int = 10) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x" * size)


def test_drops_search_and_rescoring_hits_and_their_caches(tmp_path):
    results = tmp_path / "results"
    _write(results / "search_results" / "a_prose.idXML", 100)
    _write(results / "percolator_results" / "a_per.idXML", 50)
    _write(results / "prediction_features" / "a_prose_pd.idXML", 20)
    _write(results / "prediction_features" / "summary.json")
    _write(results / "psm_filter" / "a_filter.idXML")
    _write(results / "quant_results" / "q.csv")
    cache = results / "insight_cache"
    for name in ("table_a_prose", "heatmap_a_prose", "seqview_a_per", "lineplot_a_per", "table_a_per__aa"):
        _write(cache / name / "data.parquet", 5)
    # Caches of the final PSMs stay, also those whose stem starts the same way
    for name in ("table_a_filter", "table_a_per_filter", "heatmap_a_filter"):
        _write(cache / name / "data.parquet", 5)

    files, freed = drop_intermediates(results)

    assert files == 3
    assert freed == 100 + 50 + 20 + 5 * 5
    assert not (results / "search_results").exists()
    assert not (results / "percolator_results").exists()
    assert sorted(p.name for p in (results / "prediction_features").iterdir()) == ["summary.json"]
    assert (results / "psm_filter" / "a_filter.idXML").exists()
    assert (results / "quant_results" / "q.csv").exists()
    assert sorted(p.name for p in cache.iterdir()) == ["heatmap_a_filter", "table_a_filter", "table_a_per_filter"]


def test_a_run_without_those_folders_loses_nothing(tmp_path):
    _write(tmp_path / "psm_filter" / "a_filter.idXML")
    assert drop_intermediates(tmp_path) == (0, 0)
    assert (tmp_path / "psm_filter" / "a_filter.idXML").exists()


def test_the_setting_and_sizes():
    assert keeps_final_only({"keep-final-only": True})
    assert not keeps_final_only({}) and not keeps_final_only(None)
    assert size_label(12_000) == "12 kB"
    assert size_label(340_000_000) == "340 MB"
    assert size_label(1_200_000_000) == "1.2 GB"
    assert size_label(12) == "12 B"


class _Log:
    def __init__(self):
        self.lines = []

    def log(self, message, level=0):
        self.lines.append(message)


def test_the_workflow_deletes_them_only_when_asked(tmp_path):
    from src.WorkflowTest import WorkflowTest

    _write(tmp_path / "search_results" / "a_prose.idXML", 100)
    wf = object.__new__(WorkflowTest)
    wf.logger = _Log()
    wf.results_dir = tmp_path
    wf.params = {}
    wf.keep_final_results_only()
    assert (tmp_path / "search_results").exists() and wf.logger.lines == []

    wf.params = {"keep-final-only": True}
    wf.keep_final_results_only()
    assert not (tmp_path / "search_results").exists()
    assert "deleted 1 search and rescoring result files" in wf.logger.lines[0]
