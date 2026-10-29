"""
Checks the CometAdapter settings in presets.json against what CometAdapter accepts.

CometAdapter (OpenMS release/3.5.0, src/topp/CometAdapter.cpp) refuses to run
without -force when the fragment bin size (2 x fragment_mass_tolerance) or the
bin offset do not fit the instrument type, and only accepts listed strings for
the unit and isotope options. A preset that breaks these rules fails at the
search step, so it is caught here instead.
"""
import json
from pathlib import Path

import pytest

PRESETS = json.loads(Path(__file__).parent.parent.joinpath("presets.json").read_text())["topp-workflow"]
COMET_INSTANCES = ("CometAdapter", "CometAdapter-TMT")

# Valid strings from CometAdapter.cpp
INSTRUMENTS = {"low_res", "high_res"}
PRECURSOR_UNITS = {"amu", "ppm", "Da"}
ISOTOPE_ERRORS = {"off", "0/1", "0/1/2", "0/1/2/3", "-8/-4/0/4/8", "-1/0/1/2/3"}

COMET_PRESETS = [
    (name, instance, preset[instance])
    for name, preset in PRESETS.items()
    for instance in COMET_INSTANCES
    if instance in preset
]


@pytest.mark.parametrize("name,instance,params", COMET_PRESETS, ids=lambda v: v if isinstance(v, str) else "")
def test_comet_preset_values_are_valid(name, instance, params):
    if "instrument" in params:
        assert params["instrument"] in INSTRUMENTS
    if "precursor_error_units" in params:
        assert params["precursor_error_units"] in PRECURSOR_UNITS
    if "isotope_error" in params:
        assert params["isotope_error"] in ISOTOPE_ERRORS
    if "fragment_error_units" in params:
        assert params["fragment_error_units"] == "Da"


@pytest.mark.parametrize(
    "name,instance,params",
    [p for p in COMET_PRESETS if "instrument" in p[2]],
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_instrument_presets_pass_comet_bin_checks(name, instance, params):
    """An instrument preset sets the fragment bins CometAdapter accepts for that instrument."""
    bin_size = 2 * params["fragment_mass_tolerance"]
    offset = params["fragment_bin_offset"]
    if params["instrument"] == "low_res":
        assert bin_size >= 0.8 and offset > 0.2
    else:
        assert bin_size <= 0.1 and offset <= 0.1


def test_instrument_presets_set_both_search_modes():
    """Instrument presets apply to the label-free and the TMT search alike."""
    for name, preset in PRESETS.items():
        instruments = [preset[i].get("instrument") for i in COMET_INSTANCES if i in preset]
        if any(instruments):
            assert all(i in preset for i in COMET_INSTANCES), name
            assert preset["CometAdapter"] == preset["CometAdapter-TMT"], name
