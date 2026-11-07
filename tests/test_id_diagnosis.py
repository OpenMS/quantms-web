"""Tests for the per-file identification diagnosis."""
import pandas as pd

# Imported inside the tests, as in test_overview.py: importing at collection
# pulls in the real streamlit before tests that mock it.


def _row(name, ms2=10000, targets=1500, error=1.0, percolator=True, profile=0.0,
         contaminants=2.0, edge=1.0):
    return {
        "File": name, "MS2 spectra": ms2, "Target PSMs": targets, "Decoy PSMs": 10,
        "ID rate (%)": round(100 * targets / ms2, 1) if ms2 else None, "Peptides": targets,
        "Mass error (ppm)": error, "Percolator": percolator, "Search (s)": 10,
        "Profile MS2 (%)": profile, "Contaminant PSMs (%)": contaminants, "At window edge (%)": edge,
    }


PARAMS = {"ProSE": {"Search:precursor:mass_tolerance_lower": 10.0,
                    "Search:precursor:mass_tolerance_upper": 10.0}}


def test_good_files_are_not_flagged():
    from src.common.id_diagnosis import diagnose_identifications

    per_file = pd.DataFrame([_row("a"), _row("b"), _row("c")])
    assert diagnose_identifications(per_file, PARAMS).empty


def test_causes_ranked_and_capped():
    from src.common.id_diagnosis import MAX_CAUSES, diagnose_identifications

    bad = _row("bad", targets=100, error=7.0, percolator=False, profile=80.0, contaminants=50.0, edge=30.0)
    flagged = diagnose_identifications(pd.DataFrame([_row("a"), bad]), PARAMS)
    assert list(flagged["File"]) == ["bad"]
    causes = flagged.iloc[0]["Likely causes"]
    assert causes.startswith("80 % of its MS2 spectra are profile data")
    assert "+7 ppm" in causes
    assert causes.count(".") >= MAX_CAUSES and "Percolator" not in causes


def test_outlier_against_the_run_median():
    from src.common.id_diagnosis import diagnose_identifications

    # 30 % is above the expected rate, but far below the others' 60 %
    rows = [_row(n, targets=6000) for n in "abc"] + [_row("d", targets=1600)]
    flagged = diagnose_identifications(pd.DataFrame(rows), PARAMS)
    assert list(flagged["File"]) == ["d"]
    assert "median 60 %" in flagged.iloc[0]["Likely causes"]


def test_no_ms2_and_unsearched_files():
    from src.common.id_diagnosis import diagnose_identifications

    empty = _row("empty", ms2=0, targets=0)
    unsearched = _row("unsearched") | {"ID rate (%)": None, "Target PSMs": None, "MS2 spectra": None}
    flagged = diagnose_identifications(pd.DataFrame([_row("a"), empty, unsearched]), PARAMS)
    assert list(flagged["File"]) == ["empty"]
    assert flagged.iloc[0]["Likely causes"] == "The file has no MS2 spectra to search."


def test_expected_rate_follows_the_search():
    from src.common.id_diagnosis import diagnose_identifications, search_kind

    assert search_kind(PARAMS) == "tryptic"
    iontrap = {"ProSE": {"Search:fragment:mass_tolerance_unit": "Da"}}
    nonspecific = {"ProSE": {"Search:enzyme": "unspecific cleavage"}}
    assert search_kind(iontrap) == "ion trap"
    assert search_kind(nonspecific) == "nonspecific"
    # 5 % is low for a tryptic search but fine without an enzyme
    per_file = pd.DataFrame([_row("a", targets=500)])
    assert len(diagnose_identifications(per_file, PARAMS)) == 1
    assert diagnose_identifications(per_file, nonspecific).empty
    assert diagnose_identifications(per_file, nonspecific | {}).empty


def test_unexplained_low_rate_points_at_the_settings():
    from src.common.id_diagnosis import diagnose_identifications

    flagged = diagnose_identifications(pd.DataFrame([_row("a", targets=800)]), PARAMS)
    assert "database, enzyme and modifications" in flagged.iloc[0]["Likely causes"]
