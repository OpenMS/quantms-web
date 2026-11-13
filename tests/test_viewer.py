"""The Viewer's catalog, its panels and the per-experiment scope."""
import ast
from pathlib import Path

from streamlit.testing.v1 import AppTest

from src.layout_manager import validate_layout
from src.viewer.catalog import CATALOG, COMPARISON_STEPS, DEFAULT_LAYOUT, LABELS, MODULES, PRESETS

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
    used = set(MODULES.values()) | {m for m, _ in COMPARISON_STEPS}
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
    assert at.markdown[0].value == "A:k_A:a | B:k_B_e1:b | A:k_A:a | B:k_B_e1:b | after:None:k:None"


def test_presets_are_valid_layouts():
    assert DEFAULT_LAYOUT == PRESETS[0]["layout"]
    assert len({p["name"] for p in PRESETS}) == len({p["label"] for p in PRESETS}) == len(PRESETS)
    for preset in PRESETS:
        assert validate_layout([preset["layout"]], CATALOG) == "", preset["name"]
        assert preset.get("help")
    assert [p["label"] for p in PRESETS] == ["Quantification", "Identification", "Pathways", "Quality control"]
    # Components no ready-made view shows: offered for views of the user's own
    shown = {n for p in PRESETS for row in p["layout"] for n in row}
    assert {c["name"] for c in CATALOG} - shown == {
        "clustered_heatmap", "psms", "psms_rescored", "psms_search", "protein_abundance",
        "run_summary", "id_quality", "protein_explorer",
    }


def test_each_summary_is_in_its_own_view_only():
    views = {p["name"]: {n for row in p["layout"] for n in row} for p in PRESETS}
    assert [v for v, names in views.items() if "dda_summary" in names] == ["identification"]
    assert [v for v, names in views.items() if "quant_summary" in names] == ["quantification"]


def test_every_view_and_component_is_explained():
    from src.viewer.teaching import COMPONENTS, VIEWS, tooltip

    assert set(VIEWS) == {p["name"] for p in PRESETS}
    assert set(COMPONENTS) == set(LABELS)
    assert "understand_results#protein-table" in tooltip("protein_table")


def test_default_view_follows_the_groups():
    from src.viewer.catalog import default_view

    assert default_view({"identification", "quantification", "groups"}) == "quantification"
    assert default_view({"identification", "quantification"}) == "identification"
    assert default_view(None) == "identification"


def test_earlier_views_are_moved_to_the_current_ones(tmp_path):
    import json

    from src.viewer.catalog import migrate_views

    path = tmp_path / "viewer.json"
    path.write_text(json.dumps({
        "current": "differential", "layout": [[["group_comparison"]]],
        "views": [
            {"name": "overview", "in_bar": True}, {"name": "identification", "in_bar": True},
            {"name": "custom-1", "label": "Mine", "rows": [["volcano"]]},
            {"name": "differential", "in_bar": True}, {"name": "pathways", "in_bar": False},
        ],
    }))
    migrate_views(path)
    saved = json.loads(path.read_text())
    assert [v["name"] for v in saved["views"]] == ["identification", "pathways", "custom-1"]
    assert "current" not in saved and "layout" not in saved
    before = path.read_text()
    migrate_views(path)  # once only
    assert path.read_text() == before


def test_abundance_profile_table():
    import pandas as pd

    from src.viewer.panels.abundance_profile import profile_table

    pivot = pd.DataFrame({"ProteinName": ["sp|P1|A_HUMAN", "sp|P2|B_HUMAN"], "s1": [4.0, 0.0], "s2": [8.0, 2.0]})
    table = profile_table(pivot, "ProteinName", {"s1": "ctrl", "s2": "drug"})
    p1 = table[table["accession"] == "sp|P1|A_HUMAN"]
    assert p1["x"].tolist() == [1.0, 2.0] and p1["y"].tolist() == [2.0, 3.0]
    assert p1["group"].tolist() == ["ctrl", "drug"]
    # An abundance of 0 is not measured
    assert table[table["accession"] == "sp|P2|B_HUMAN"]["sample"].tolist() == ["s2"]
