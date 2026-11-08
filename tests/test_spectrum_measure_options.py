"""The PSM spectrum's measure tool knows the PSM's modified residues and the search's tolerance."""

from src.common.spectrum_view import measure_options, modified_residues


def test_modified_residues_of_a_psm():
    masses = modified_residues.__wrapped__(".(Acetyl)PEPM(Oxidation)C(Carbamidomethyl)K")
    # The N-terminal acetyl has no residue of its own; the two modified residues do
    assert set(masses) == {"M(Oxidation)", "C(Carbamidomethyl)"}
    assert abs(masses["M(Oxidation)"] - 147.0354) < 1e-4
    assert abs(masses["C(Carbamidomethyl)"] - 160.03065) < 1e-4


def test_unmodified_or_unreadable_peptides_give_nothing():
    assert modified_residues.__wrapped__("PEPTIDEK") == {}
    assert modified_residues.__wrapped__("PEP(NotAModification)K") == {}


def test_tolerance_in_da_reaches_the_plot_ppm_does_not():
    assert measure_options("PEPTIDEK", 0.02, False) == {"measure_residues": {}, "measure_tolerance": 0.02}
    assert measure_options("PEPTIDEK", 20.0, True) == {"measure_residues": {}}
