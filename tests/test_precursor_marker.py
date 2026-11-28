"""The spectra mark the precursor with the visualization library's marker.

The modules are imported inside the tests, not at collection: importing
spectrum_view then would load src.workflow before the tests that mock
streamlit for it are collected.
"""
from types import SimpleNamespace


def test_marker_comes_from_the_library(monkeypatch):
    import src.viz.identification.annotated_spectrum as library
    from src.common import spectrum_view

    calls = []
    monkeypatch.setattr(library, "precursor_marker", lambda cls, mz, charge: calls.append((cls, mz, charge)) or {"x_markers": []}, raising=False)
    assert spectrum_view.precursor_marker(SimpleNamespace, 512.3, 2) == {"x_markers": []}
    assert calls == [(SimpleNamespace, 512.3, 2)]


def test_no_marker_with_a_library_without_it(monkeypatch):
    import src.viz.identification.annotated_spectrum as library
    from src.common import spectrum_view

    monkeypatch.delattr(library, "precursor_marker", raising=False)
    assert spectrum_view.precursor_marker(SimpleNamespace, 512.3, 2) == {}
