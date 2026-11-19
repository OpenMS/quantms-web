"""Keep only the final results of a run: drop the search and rescoring hits.

A run keeps every step's idXML files. The raw search hits and the rescored
hits of a large study take much of its disk, and only the views that show
them (Raw search hits, Rescored hits, score distributions, identification
quality) need them. With **Keep only final results** the workflow deletes
them, and the viewer caches made from them, once the run has finished; the
filtered PSMs, the proteins and the abundances stay. Streamlit-free, so the
workflow process can call it.
"""
from __future__ import annotations

import shutil
from pathlib import Path

# The workflow setting (Describe experiment, False discovery rate tab)
FINAL_ONLY_KEY = "keep-final-only"
# Folders of results/ deleted whole: raw search hits and rescored hits
DROPPED_FOLDERS = ("search_results", "percolator_results")
# Folders of results/ whose idXML files are deleted; the rest (summaries) stay
DROPPED_IDXML = ("prediction_features",)
CACHE_DIR = "insight_cache"
# The caches src/common/psm_views.py makes per idXML stem
CACHE_PREFIXES = ("table", "heatmap", "seqview", "lineplot")


def keeps_final_only(params: dict) -> bool:
    """Whether a run's settings ask to keep only its final results."""
    return bool((params or {}).get(FINAL_ONLY_KEY, False))


def _size(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def _drop_caches(cache_dir: Path, stem: str) -> int:
    """Delete the viewer caches of one idXML; the bytes freed."""
    if not cache_dir.is_dir():
        return 0
    targets = [cache_dir / f"{prefix}_{stem}" for prefix in CACHE_PREFIXES]
    targets += list(cache_dir.glob(f"table_{stem}__*"))
    freed = 0
    for target in targets:
        if target.exists():
            freed += _size(target)
            shutil.rmtree(target) if target.is_dir() else target.unlink()
    return freed


def drop_intermediates(results_dir) -> tuple[int, int]:
    """Delete the search and rescoring hits of a run and their viewer caches.

    Returns:
        ``(files, bytes)``: the idXML files deleted and the disk freed.
    """
    results = Path(results_dir)
    cache_dir = results / CACHE_DIR
    files = freed = 0
    for name in DROPPED_FOLDERS + DROPPED_IDXML:
        folder = results / name
        if not folder.is_dir():
            continue
        idxmls = sorted(folder.glob("*.idXML"))
        for idxml in idxmls:
            freed += _drop_caches(cache_dir, idxml.stem)
        files += len(idxmls)
        if name in DROPPED_FOLDERS:
            freed += _size(folder)
            shutil.rmtree(folder)
        else:
            for idxml in idxmls:
                freed += idxml.stat().st_size
                idxml.unlink()
    return files, freed


def size_label(n_bytes: int) -> str:
    """'1.2 GB', '340 MB', '12 kB'."""
    for unit, scale in (("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if n_bytes >= scale:
            return f"{n_bytes / scale:.1f} {unit}" if n_bytes < 10 * scale else f"{n_bytes / scale:.0f} {unit}"
    return f"{n_bytes} B"
