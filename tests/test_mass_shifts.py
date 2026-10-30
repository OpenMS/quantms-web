"""Tests for the precursor mass shifts and their Unimod interpretation (Mass Shifts page)."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyopenms as poms

# src.common.mass_shifts is imported inside the tests, not at collection: it pulls in
# src.workflow.ParameterManager with the real streamlit, and test_parameter_defaults
# needs to import it with a mock.

PROTON = 1.007276466812


def _idxml(path, hits):
    """hits: (sequence, charge, mass shift in Da, "target" or "decoy")."""
    prot = poms.ProteinIdentification()
    prot.setIdentifier("run")
    hit_p = poms.ProteinHit()
    hit_p.setAccession("P1")
    prot.setHits([hit_p])
    peptides = []
    for seq, charge, shift, td in hits:
        aa = poms.AASequence.fromString(seq)
        pep = poms.PeptideIdentification()
        pep.setIdentifier("run")
        pep.setMZ((aa.getMonoWeight() + shift + charge * PROTON) / charge)
        hit = poms.PeptideHit()
        hit.setSequence(aa)
        hit.setCharge(charge)
        hit.setMetaValue("target_decoy", td)
        evidence = poms.PeptideEvidence()
        evidence.setProteinAccession("P1")
        hit.setPeptideEvidences([evidence])
        pep.setHits([hit])
        peptides.append(pep)
    if hasattr(poms, "PeptideIdentificationList"):
        container = poms.PeptideIdentificationList()
        for pep in peptides:
            container.push_back(pep)
        peptides = container
    path.parent.mkdir(parents=True, exist_ok=True)
    poms.IdXMLFile().store(str(path), [prot], peptides)


def test_psm_mass_shifts_reads_targets_with_their_shift(tmp_path):
    from src.common.mass_shifts import psm_mass_shifts

    _idxml(tmp_path / "runA_filter.idXML", [
        ("PEPTIDEK", 2, 0.0, "target"),
        ("PEPTIDEK", 3, 15.9949, "target"),
        ("EDITPEPK", 2, 0.0, "decoy"),
    ])
    psms = psm_mass_shifts([tmp_path / "runA_filter.idXML"], ["runA.mzML"])
    assert len(psms) == 2
    assert set(psms["file"]) == {"runA.mzML"}
    assert np.allclose(sorted(psms["delta"]), [0.0, 15.9949], atol=1e-4)


def test_unimod_matches_put_common_modifications_first():
    from src.common.mass_shifts import unimod_matches

    matches = unimod_matches(79.9663)
    assert matches[0]["name"] == "Phospho"
    assert set("STY") <= set(matches[0]["residues"])
    assert abs(matches[0]["deviation"]) < 1e-3
    assert not any(m["name"].startswith("MOD:") for m in matches)


def test_summary_finds_and_interprets_shifts():
    from src.common.mass_shifts import mass_shift_summary

    rng = np.random.default_rng(1)
    deltas = np.concatenate([
        rng.normal(0.0, 0.002, 40),
        rng.normal(15.9949, 0.002, 12),
        rng.normal(1.00335, 0.002, 6),
        rng.normal(-17.0265, 0.002, 4),
        [123.4, 250.7],  # lone shifts: below the minimum PSM count
    ])
    psms = pd.DataFrame({
        "file": "a.mzML", "sequence": "PEPTIDEK", "charge": 2, "mz": 500.0,
        "rt": 1.0, "score": 0.01, "delta": deltas,
    })

    shifts, histogram, linked = mass_shift_summary(psms)

    # Most PSMs first, ids in that order
    assert shifts["shift_id"].tolist() == [0, 1, 2, 3]
    assert shifts["PSMs"].tolist() == [40, 12, 6, 4]
    assert shifts["Interpretation"].tolist() == [
        "Unmodified", "Oxidation", "13C isotope peak", "Ammonia-loss",
    ]
    assert abs(shifts.loc[1, "Deviation (mDa)"]) < 3
    assert "M" in shifts.loc[1, "Residues"]

    # Every PSM is linked to its shift; the lone ones to none
    assert (linked["shift_id"] == -1).sum() == 2
    assert linked.groupby("shift_id").size().to_dict() == {-1: 2, 0: 40, 1: 12, 2: 6, 3: 4}

    # One labelled peak per shift in the histogram
    peaks = histogram[histogram["peak"]]
    assert sorted(peaks["shift_id"]) == [0, 1, 2, 3]
    assert peaks.set_index("shift_id").loc[1, "label"].endswith("Oxidation")
    assert histogram["psms"].sum() == len(psms)


def test_deamidation_and_isotope_pick_stay_apart():
    from src.common.mass_shifts import mass_shift_summary

    rng = np.random.default_rng(2)
    # Spread enough that the two peaks touch
    deltas = np.concatenate([rng.normal(0.984, 0.003, 30), rng.normal(1.00335, 0.003, 30)])
    psms = pd.DataFrame({
        "file": "a.mzML", "sequence": "PEPTIDEK", "charge": 2, "mz": 500.0,
        "rt": 1.0, "score": 0.01, "delta": deltas,
    })
    shifts, _, _ = mass_shift_summary(psms)
    assert sorted(shifts["Interpretation"]) == ["13C isotope peak", "Deamidated"]
    assert all(abs(n - 30) <= 2 for n in shifts["PSMs"])


def test_one_noisy_peak_is_not_split():
    from src.common.mass_shifts import mass_shift_summary

    deltas = np.random.default_rng(3).normal(15.9949, 0.003, 500)
    psms = pd.DataFrame({
        "file": "a.mzML", "sequence": "PEPTIDEK", "charge": 2, "mz": 500.0,
        "rt": 1.0, "score": 0.01, "delta": deltas,
    })
    shifts, _, _ = mass_shift_summary(psms)
    assert shifts["PSMs"].tolist() == [500]


def test_unexplained_shift_has_no_interpretation():
    from src.common.mass_shifts import mass_shift_summary

    psms = pd.DataFrame({
        "file": "a.mzML", "sequence": "PEPTIDEK", "charge": 2, "mz": 500.0,
        "rt": 1.0, "score": 0.01, "delta": [333.333] * 5,
    })
    shifts, _, _ = mass_shift_summary(psms)
    assert shifts.loc[0, "Interpretation"] == ""
    assert pd.isna(shifts.loc[0, "Unimod Δ (Da)"])


def test_open_search_preset_is_valid_comet_parameters():
    presets = json.loads(Path("presets.json").read_text())["topp-workflow"]
    for name, expected_units in (("Open Search", "Da"), ("Closed Search", "ppm")):
        for tool in ("CometAdapter", "CometAdapter-TMT"):
            params = presets[name][tool]
            assert params["precursor_error_units"] == expected_units
            # CometAdapter's valid strings
            assert params["precursor_error_units"] in ("amu", "ppm", "Da")
            assert params["isotope_error"] in ("off", "0/1", "0/1/2", "0/1/2/3", "-8/-4/0/4/8", "-1/0/1/2/3")
    assert presets["Open Search"]["CometAdapter"]["precursor_mass_tolerance"] == 500.0
