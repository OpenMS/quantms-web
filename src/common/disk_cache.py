"""Results the Viewer reads, kept on disk with the run they were computed from.

The workflow computes them once at the end of a run (``src/viewer/precompute.py``),
so the Viewer opens on what is already there instead of parsing idXMLs or
comparing groups on its first visit. A result is kept under ``<run>/viewer-cache/``
(the nearest folder holding ``run.json``, else the workflow folder), keyed on the
function, its arguments and the files it reads (``file_id``: where the file sits
in its run and what it holds): a new run or a rewritten file is computed again.

The keys name no absolute path and no modification time, so a run copied
elsewhere keeps its results: a demo workspace ships its ``viewer-cache/``, and
every workspace set up from it opens the Viewer on it. Streamlit-free.
"""
from __future__ import annotations

import functools
import hashlib
import os
import pickle
from pathlib import Path
from typing import Any, Callable

FOLDER = "viewer-cache"

# Bytes read from the start, middle and end of a file for its fingerprint
_SAMPLE = 64 * 1024
_FINGERPRINTS: dict = {}


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


def _fingerprint(path: Path, stat: os.stat_result) -> str:
    """A hash of a file's size and of samples from its start, middle and end.

    Kept per process for the file's modification time, so a rerun of the page
    does not read the file again.
    """
    memo = (str(path.absolute()), stat.st_mtime_ns, stat.st_size)
    if memo not in _FINGERPRINTS:
        digest = hashlib.sha256(str(stat.st_size).encode())
        with open(path, "rb") as handle:
            for offset in sorted({0, max(0, stat.st_size // 2 - _SAMPLE // 2), max(0, stat.st_size - _SAMPLE)}):
                handle.seek(offset)
                digest.update(handle.read(_SAMPLE))
        _FINGERPRINTS[memo] = digest.hexdigest()[:24]
    return _FINGERPRINTS[memo]


def file_id(path, root=None) -> tuple:
    """A file as a cache key: its path within ``root`` (else its name), size and fingerprint.

    The same in every copy of a run, linked or copied, whatever its modification time.
    """
    path = Path(path)
    stat = path.stat()
    try:
        name = path.absolute().relative_to(Path(root).absolute()).as_posix() if root else path.name
    except ValueError:
        name = path.name
    return name, stat.st_size, _fingerprint(path, stat)


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
    """Keep ``func(path, *args)`` on disk, keyed on the file (``file_id``) and the arguments."""

    @functools.wraps(func)
    def wrapper(path, *args, **kwargs):
        try:
            root = cache_root(Path(path))
            ident = file_id(path, root)
        except (OSError, TypeError):
            return func(path, *args, **kwargs)
        key = (func.__module__, func.__qualname__, ident, args, sorted(kwargs.items()))
        return remember(root, func.__name__, key, lambda: func(path, *args, **kwargs))

    wrapper.uncached = func
    return wrapper
