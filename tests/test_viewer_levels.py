"""The protein, peptide and PSM tables take their place from the library."""
from types import SimpleNamespace

from src.viewer.protein_data import _place


def test_place_comes_from_the_library_display():
    seen = {}

    def display(cls, *args):
        seen["cls"], seen["args"] = cls, args
        return {"level": 1}

    assert _place(SimpleNamespace(display=display), "protein") == {"level": 1}
    assert seen["cls"].__name__ == "Table"
    assert seen["args"] == ("protein",)


def test_place_is_empty_with_a_library_without_it():
    assert _place(SimpleNamespace()) == {}
