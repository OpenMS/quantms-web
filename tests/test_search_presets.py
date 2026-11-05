"""
Checks the ProSE settings in presets.json against what ProSE accepts.

ProSE (OpenMS 3.6.0, src/topp/ProSE.cpp) only accepts "ppm" or "Da" as tolerance
units and integer isotope errors (stored one per line, as the settings page writes
list settings), and switches to an open search when the
precursor tolerance is above 1 Da. A preset that breaks these rules fails at the
search step, or searches in a mode it does not name, so it is caught here instead.
"""
import json
from pathlib import Path

import pytest

PRESETS = json.loads(Path(__file__).parent.parent.joinpath("presets.json").read_text())["topp-workflow"]
SEARCH_INSTANCES = ("ProSE", "ProSE-TMT")
UNITS = {"ppm", "Da"}

SEARCH_PRESETS = [
    (name, instance, preset[instance])
    for name, preset in PRESETS.items()
    for instance in SEARCH_INSTANCES
    if instance in preset
]


def test_presets_set_only_search_settings():
    """Every preset names a search instance, and no tool the workflow no longer runs."""
    for name, preset in PRESETS.items():
        tools = {k for k in preset if not k.startswith("_")}
        assert tools and tools <= set(SEARCH_INSTANCES), name


@pytest.mark.parametrize("name,instance,params", SEARCH_PRESETS, ids=lambda v: v if isinstance(v, str) else "")
def test_search_preset_values_are_valid(name, instance, params):
    assert all(key.startswith("Search:") for key in params)
    for key in ("Search:precursor:mass_tolerance_unit", "Search:fragment:mass_tolerance_unit"):
        if key in params:
            assert params[key] in UNITS
    if "Search:precursor:isotopes" in params:
        # A list setting is stored one value per line, as its text field writes it
        isotopes = params["Search:precursor:isotopes"]
        assert isinstance(isotopes, str) and all(i.lstrip("-").isdigit() for i in isotopes.split("\n"))
    # The precursor window is symmetric
    assert params.get("Search:precursor:mass_tolerance_lower") == params.get("Search:precursor:mass_tolerance_upper")


@pytest.mark.parametrize("name,instance,params", SEARCH_PRESETS, ids=lambda v: v if isinstance(v, str) else "")
def test_open_search_only_where_named(name, instance, params):
    """ProSE searches openly above 1 Da (or 1000 ppm); only the Open Search preset does that."""
    tolerance = params.get("Search:precursor:mass_tolerance_upper")
    if tolerance is None:
        return
    unit = params.get("Search:precursor:mass_tolerance_unit", "ppm")
    is_open = tolerance > (1.0 if unit == "Da" else 1000.0)
    assert is_open == (name == "Open Search"), name


def test_instrument_presets_set_both_search_modes():
    """Instrument presets apply to the label-free and the TMT search alike."""
    for name, preset in PRESETS.items():
        if "Search:fragment:mass_tolerance" in preset.get("ProSE", {}):
            assert preset.get("ProSE") == preset.get("ProSE-TMT"), name
