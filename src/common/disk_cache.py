"""Results the Viewer reads, kept on disk with the run they were computed from.

The workflow computes them once at the end of a run (``src/viewer/precompute.py``),
so the Viewer opens on what is already there instead of parsing idXMLs or
comparing groups on its first visit. A result is kept under ``<run>/viewer-cache/``
(the nearest folder holding ``run.json``, else the workflow folder), keyed on the
function, its arguments and the size and modification time of the file it reads:
a new run or a rewritten file is computed again. Streamlit-free.
"""
from __future__ import annotations

import functools
import hashlib
import os
import pickle
from pathlib import Path
from typing import Any, Callable

FOLDER = "viewer-cache"


def cache_root(path: Path) -> Path:
    """The run folder ``path`` belongs to (``run.json``), else its workflow folder, else its folder."""
    path = Path(path)
    for parent in path.parents:
        if (parent / "run.json").exists():
            return parent
    for parent in path.parents:
        if (parent / "params.json").exists():
            return parent
    return path.parent


def _key(*parts: Any) -> str:
    return hashlib.sha256(repr(parts).encode()).hexdigest()[:24]


def remember(root: Path, name: str, key: Any, compute: Callable[[], Any]) -> Any:
    """``compute()``, kept as ``<root>/viewer-cache/<name>-<key>.pkl``.

    A file that cannot be read is computed again; one that cannot be written
    (a read-only folder) is not kept.
    """
    target = Path(root) / FOLDER / f"{name}-{_key(key)}.pkl"
    if target.is_file():
        try:
            with open(target, "rb") as handle:
                return pickle.load(handle)
        except Exception:  # pylint: disable=broad-except
            pass
    value = compute()
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        partial = target.with_suffix(f".{os.getpid()}.part")
        with open(partial, "wb") as handle:
            pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(partial, target)
    except OSError:
        pass
    return value


def per_file(func: Callable) -> Callable:
    """Keep ``func(path, *args)`` on disk, keyed on the file's size and modification time."""

    @functools.wraps(func)
    def wrapper(path, *args, **kwargs):
        try:
            stat = Path(path).stat()
        except (OSError, TypeError):
            return func(path, *args, **kwargs)
        key = (func.__module__, func.__qualname__, str(Path(path).resolve()), stat.st_mtime_ns, stat.st_size,
               args, sorted(kwargs.items()))
        return remember(cache_root(Path(path)), func.__name__, key, lambda: func(path, *args, **kwargs))

    wrapper.uncached = func
    return wrapper
