"""The Viewer's catalog, its panels and the per-experiment scope."""
import ast
from pathlib import Path

from streamlit.testing.v1 import AppTest

from src.viewer.catalog import CATALOG, COMPARISON_STEPS, DEFAULT_LAYOUT, LABELS, MODULES

PANELS = Path("src/viewer/panels")


def _has_render(module: str) -> bool:
    tree = ast.parse((PANELS / f"{module}.py").read_text(encoding="utf-8"))
    return any(isinstance(n, ast.FunctionDef) and n.name == "render" for n in tree.body)


def test_every_component_has_a_panel():
    for name, module in MODULES.items():
        assert (PANELS / f"{module}.py").is_file(), name
        assert _has_render(module), module


def test_every_comparison_step_has_a_panel():
    for module, _ in COMPARISON_STEPS:
        assert _has_render(module), module


def test_names_are_unique_and_labelled():
    names = [c["name"] for c in CATALOG]
    assert len(names) == len(set(names))
    assert all(LABELS[n] for n in names)


def test_default_layout_uses_catalog_components():
    for row in DEFAULT_LAYOUT:
        assert 1 <= len(row) <= 3
        assert all(name in MODULES for name in row)


def test_every_panel_is_used():
    used = set(MODULES.values()) | {m for m, _ in COMPARISON_STEPS} | {"download", "library"}
    on_disk = {p.stem for p in PANELS.glob("*.py") if p.stem != "__init__"}
    assert on_disk == used


_SCOPE_APP = """
import streamlit as st
from src.viewer.scope import experiment, scoped, shown_run

out = []
st.session_state.setdefault("_viewer_chains", {})
for run, index, value in [("A", 0, "a"), ("B", 1, "b"), ("A", 0, None), ("B", 1, None)]:
    with experiment(run, index):
        if value is not None:
            st.session_state["filtered_df"] = value
        out.append(f"{shown_run()}:{scoped('k')}:{st.session_state.get('filtered_df')}")
out.append(f"after:{shown_run()}:{scoped('k')}:{st.session_state.get('filtered_df')}")
st.markdown(" | ".join(out))
"""


def test_experiments_keep_their_own_comparison():
    at = AppTest.from_string(_SCOPE_APP)
    at.run(timeout=30)
    assert not at.exception
    assert at.markdown[0].value == "A:k:a | B:k_e1:b | A:k:a | B:k_e1:b | after:None:k:None"
