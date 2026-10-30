"""Tests for the peptide display and export forms (labels, ProForma)."""
import pytest

from src.common.peptide_format import format_sequence, to_labels, to_proforma

MODIFIED = ".(TMT6plex)PEPM(Oxidation)TIDEC(Carbamidomethyl)S(Phospho)K(TMT6plex)"


def test_proforma_uses_unimod_names_and_terminal_dashes():
    assert to_proforma(MODIFIED) == "[TMT6plex]-PEPM[Oxidation]TIDEC[Carbamidomethyl]S[Phospho]K[TMT6plex]"
    assert to_proforma("PEPTIDEK.(Amidated)") == "PEPTIDEK-[Amidated]"
    assert to_proforma("PEPTIDEK") == "PEPTIDEK"


def test_labels_shorten_known_modifications():
    assert to_labels(MODIFIED) == "tmt-PEPM(ox)TIDEC(cam)S(ph)K(tmt)"
    # A mass-only modification resolves to its Unimod name
    assert to_labels("PEPN[+0.984]K") == "PEPN(de)K"


def test_mass_shift_without_a_name():
    # No Unimod entry at this mass, so ProForma carries the signed mass shift
    proforma = to_proforma("PEPK[+123.4567]")
    assert proforma.startswith("PEPK[+") and proforma.endswith("]")
    assert abs(float(proforma[5:-1]) - 123.4567) < 1e-3


@pytest.mark.parametrize("mode, expected", [
    ("Labels", "PEPM(ox)K"),
    ("ProForma", "PEPM[Oxidation]K"),
    ("OpenMS", "PEPM(Oxidation)K"),
])
def test_format_sequence_modes(mode, expected):
    assert format_sequence("PEPM(Oxidation)K", mode) == expected


def test_unreadable_strings_pass_through():
    assert to_proforma("not a peptide!") == "not a peptide!"
    assert to_labels("not a peptide!") == "not a peptide!"
