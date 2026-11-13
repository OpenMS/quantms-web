"""The AlphaPeptDeep predictions run outside the app server (torch must not load there)."""
import io
import json
import sys
import textwrap

import pandas as pd
import pytest

from src.common import prediction_worker
from src.common.prediction_worker import PredictionError, PredictionWorker


def test_peptdeep_available_does_not_import_it(monkeypatch):
    import importlib.util

    import src.common.ms2_prediction as ms2

    monkeypatch.delitem(sys.modules, "peptdeep", raising=False)
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: object() if name == "peptdeep" else None)
    assert ms2.peptdeep_available()
    assert "peptdeep" not in sys.modules
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None)
    assert not ms2.peptdeep_available()


def test_serve_answers_each_request_and_survives_a_failure(monkeypatch):
    import src.common.ms2_prediction as ms2

    loads = []

    def predict(model, peptide, charge, instrument):
        if peptide == "BAD":
            raise RuntimeError("no such residue")
        if peptide == "SHORT":
            return None
        return pd.DataFrame({"mz": [147.1], "intensity": [1.0], "ion": ["y1¹⁺"], "extra": [0]})

    monkeypatch.setattr(ms2, "load_model", lambda: loads.append(1) or "model")
    monkeypatch.setattr(ms2, "predict_fragments", predict)
    requests = io.StringIO("".join(
        json.dumps({"peptide": p, "charge": 2, "instrument": "Lumos"}) + "\n" for p in ["PEPTIDEK", "BAD", "SHORT"]
    ))
    answers = io.StringIO()
    prediction_worker.serve(requests, answers)
    lines = [json.loads(line) for line in answers.getvalue().splitlines()]
    assert lines[0] == {"rows": [{"mz": 147.1, "intensity": 1.0, "ion": "y1¹⁺"}]}
    assert "no such residue" in lines[1]["error"]
    assert lines[2] == {"rows": None}
    assert loads == [1]  # the model is loaded once


FAKE_CHILD = textwrap.dedent("""
    import json, sys
    print("noise on stderr", file=sys.stderr)
    for line in sys.stdin:
        request = json.loads(line)
        if request["peptide"] == "DIE":
            sys.exit(1)
        if request["peptide"] == "BAD":
            print(json.dumps({"error": "bad peptide"}), flush=True)
            continue
        rows = [{"mz": 100.0 + request["charge"], "intensity": 1.0, "ion": "b1¹⁺"}]
        print(json.dumps({"rows": rows}), flush=True)
""")


def test_worker_talks_to_its_process_and_restarts_it(tmp_path):
    script = tmp_path / "child.py"
    script.write_text(FAKE_CHILD, encoding="utf-8")
    worker = PredictionWorker([sys.executable, str(script)])
    try:
        predicted = worker.predict("PEPTIDEK", 2, "Lumos")
        assert predicted.to_dict("records") == [{"mz": 102.0, "intensity": 1.0, "ion": "b1¹⁺"}]
        with pytest.raises(PredictionError, match="bad peptide"):
            worker.predict("BAD", 2, "Lumos")
        with pytest.raises(PredictionError, match="stopped"):
            worker.predict("DIE", 2, "Lumos")
        # the next request starts a new process
        assert worker.predict("PEPTIDEK", 3, "Lumos")["mz"].tolist() == [103.0]
    finally:
        worker.close()
