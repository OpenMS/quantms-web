"""Single run | Comparison: a run shown on its own compares its own files only."""
from types import SimpleNamespace

import streamlit as st

from src.common import results_helpers
from src.viewer import scope


def test_experiment_runs_single_run(monkeypatch):
    runs = SimpleNamespace(experiment=lambda run_id: ["r1", "r2"])
    monkeypatch.setattr(results_helpers, "get_runs", lambda workspace: runs)
    st.session_state.clear()
    with scope.experiment("r1", 0):
        assert results_helpers.experiment_runs("ws") == ["r1", "r2"]
    with scope.experiment("r2", 0, single=True):
        assert results_helpers.experiment_runs("ws") == ["r2"]
        # Another run's experiment is still the whole of it
        assert results_helpers.experiment_runs("ws", "r1") == ["r1", "r2"]
    assert not scope.shown_single()


def test_single_run_and_whole_keep_comparisons_apart():
    st.session_state.clear()
    from src.common.results_helpers import RUN_DERIVED_KEYS

    key = RUN_DERIVED_KEYS[0]
    with scope.experiment("r1", 0):
        st.session_state[key] = "whole"
    with scope.experiment("r1", 0, single=True):
        assert key not in st.session_state
        st.session_state[key] = "single"
    with scope.experiment("r1", 0):
        assert st.session_state[key] == "whole"
