import pandas as pd
import pytest
from pyopenms import AASequence, Residue

from src.common.fragment_plots import fragment_error_plot, ion_table, ion_table_html, load_psm_fragments


def _fragments():
    return pd.DataFrame({
        "sequence_id": [0, 0, 0, 0],
        "sequence": ["PEPTIDEK"] * 4,
        "precursor_charge": [2] * 4,
        "ion": ["b2", "b2-H2O", "y3", "y3"],
        "ion_type": ["b", "b", "y", "y"],
        "ion_number": [2, 2, 3, 3],
        "charge": [1, 1, 1, 1],
        "theoretical_mz": [227.1, 209.1, 391.2, 391.2],
        "observed_mz": [227.1, 209.1, 391.2, 391.21],
        "intensity": [100.0, 50.0, 200.0, 20.0],
        "mass_error_da": [0.001, -0.002, 0.0005, 0.01],
        "mass_error_ppm": [4.0, -9.0, 1.2, 25.0],
        "peak_id": [1, 2, 3, 4],
    })


def test_ion_table_values_and_matches():
    values, matched = ion_table("PEPTIDEK", 2, _fragments())
    assert list(values.columns) == ["#", "b+", "residue", "y+", "#y"]
    peptide = AASequence.fromString("PEPTIDEK")
    assert values.loc[1, "b+"] == pytest.approx(peptide.getPrefix(2).getMonoWeight(Residue.ResidueType.BIon, 1))
    # row i holds y_(n-i+1): y3 sits in row 6 (index 5)
    assert values.loc[5, "#y"] == 3
    assert values.loc[5, "y+"] == pytest.approx(peptide.getSuffix(3).getMonoWeight(Residue.ResidueType.YIon, 1))
    assert matched.loc[1, "b+"] and matched.loc[5, "y+"]
    assert matched["b+"].sum() == 1  # b2-H2O is a neutral loss, not counted
    assert pd.isna(values.loc[7, "b+"]) and pd.isna(values.loc[0, "y+"])


def test_ion_table_doubly_charged_columns():
    values, _ = ion_table("PEPTIDEK", 3, _fragments())
    assert list(values.columns) == ["#", "b+", "b++", "residue", "y+", "y++", "#y"]


def test_ion_table_html_highlights_matches():
    values, matched = ion_table("PEPTIDEK", 2, _fragments())
    out = ion_table_html(values, matched)
    assert out.count('class="hit-b"') == 1 and out.count('class="hit-y"') == 1


def test_fragment_error_plot_uses_tolerance_unit():
    fig = fragment_error_plot(_fragments(), 20, True)
    assert fig.layout.yaxis.title.text == "Mass error (ppm)"
    assert [t.name for t in fig.data] == ["b ions", "y ions"]
    fig_da = fragment_error_plot(_fragments(), 0.02, False)
    assert list(fig_da.data[0].y) == [0.001, -0.002]


def test_load_psm_fragments_filters_one_psm(tmp_path):
    tsv = tmp_path / "run_fragment_ions.tsv"
    frags = _fragments()
    pd.concat([frags, frags.assign(sequence_id=1)]).to_csv(tsv, sep="\t", index=False)
    assert len(load_psm_fragments(tsv, 1)) == 4
    assert load_psm_fragments(tsv, 7).empty
