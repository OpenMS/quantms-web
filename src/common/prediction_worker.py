"""AlphaPeptDeep predictions in a process of their own.

The Viewer mirrors a PSM's spectrum against its predicted spectrum. peptdeep
loads torch, and on Windows polars' parquet scan crashed the app server
(access violation) once torch was loaded in the same process. So the model
lives in a child process, ``python -m src.common.prediction_worker``, that
loads it once and answers one JSON line per request on stdin with one JSON line
on stdout. Streamlit-free.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import threading
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[2]


class PredictionError(RuntimeError):
    """The prediction process failed or answered with an error."""


class PredictionWorker:
    """Client of one prediction process, started on first use and again when it died."""

    def __init__(self, command: list[str] | None = None):
        self._command = command or [sys.executable, "-m", "src.common.prediction_worker"]
        self._process: subprocess.Popen | None = None
        self._lock = threading.Lock()

    def _start(self) -> subprocess.Popen:
        flags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
        return subprocess.Popen(
            self._command, cwd=_ROOT, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=None, text=True, encoding="utf-8", bufsize=1, creationflags=flags,
        )

    def predict(self, peptide: str, charge: int, instrument: str) -> pd.DataFrame | None:
        """Predicted b/y ions (columns mz, intensity, ion), None when the peptide cannot be predicted."""
        request = json.dumps({"peptide": peptide, "charge": int(charge), "instrument": instrument})
        with self._lock:
            if self._process is None or self._process.poll() is not None:
                self._process = self._start()
            try:
                self._process.stdin.write(request + "\n")
                self._process.stdin.flush()
                line = self._process.stdout.readline()
            except (BrokenPipeError, OSError) as error:
                self._process = None
                raise PredictionError(f"the prediction process stopped: {error}") from error
            if not line:
                self._process = None
                raise PredictionError("the prediction process stopped without an answer")
        answer = json.loads(line)
        if "error" in answer:
            raise PredictionError(answer["error"])
        rows = answer.get("rows")
        return None if rows is None else pd.DataFrame(rows, columns=["mz", "intensity", "ion"])

    def close(self) -> None:
        with self._lock:
            if self._process is not None and self._process.poll() is None:
                self._process.stdin.close()
                try:
                    self._process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    self._process.kill()
            self._process = None


def serve(requests=None, answers=None) -> None:
    """Answer prediction requests until stdin closes."""
    if answers is None:
        # Anything the libraries print goes to stderr, so stdout carries only answers.
        answers = os.fdopen(os.dup(sys.stdout.fileno()), "w", encoding="utf-8", buffering=1)
        os.dup2(sys.stderr.fileno(), sys.stdout.fileno())
        sys.stdout = sys.stderr
    requests = requests if requests is not None else sys.stdin
    from src.common.ms2_prediction import load_model, predict_fragments

    model = None
    for line in requests:
        if not line.strip():
            continue
        try:
            request = json.loads(line)
            if model is None:
                model = load_model()
            predicted = predict_fragments(model, request["peptide"], int(request["charge"]), request["instrument"])
            answer = {"rows": None if predicted is None else predicted[["mz", "intensity", "ion"]].to_dict("records")}
        except Exception as error:  # one failed request must not end the process
            answer = {"error": f"{type(error).__name__}: {error}"}
        answers.write(json.dumps(answer) + "\n")
        answers.flush()


if __name__ == "__main__":
    serve()
