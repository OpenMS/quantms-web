"""File name of a PSM spectrum's SVG download."""


def test_plain_peptide():
    from src.common.spectrum_view import svg_name

    assert svg_name("PEPTIDEK", 2, 1234) == "PEPTIDEK_2+_scan1234"


def test_modifications_become_safe_characters():
    from src.common.spectrum_view import svg_name

    assert svg_name("PEPM(Oxidation)TIDEK", 3, 7) == "PEPM-Oxidation-TIDEK_3+_scan7"
    assert svg_name(".(Acetyl)PEPTIDEK", 2, 1) == "Acetyl-PEPTIDEK_2+_scan1"
    assert svg_name("PEPC[160.03]K", 2, 1) == "PEPC-160.03-K_2+_scan1"


def test_nothing_unsafe_reaches_the_file_name():
    from src.common.spectrum_view import svg_name

    name = svg_name('A/B\\C:D*E?F"G<H>I|K', 2, 5)
    assert not set('/\\:*?"<>| ') & set(name)
    assert svg_name("", 2, 5) == "spectrum_2+_scan5"
