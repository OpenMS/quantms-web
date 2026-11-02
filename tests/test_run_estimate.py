"""Tests for src/common/run_estimate.py."""
from src.common.run_estimate import (
    COMET_BASE_MB,
    COMET_MB_PER_SPECTRUM,
    GB,
    estimate,
    fasta_stats,
    format_minutes,
    search_factor,
    spectrum_count,
)

REFERENCE = {
    "num_enzyme_termini": "fully", "variable_modifications": "Oxidation (M)\nAcetyl (Protein N-term)",
    "max_variable_mods_in_peptide": "3", "missed_cleavages": "2", "precursor_mass_tolerance": "20.0",
    "precursor_error_units": "ppm", "instrument": "high_res",
}


def _mzml(path, count, size=None):
    path.write_text(f'<?xml version="1.0"?><mzML><run><spectrumList count="{count}" defaultDataProcessingRef="x">')
    if size:
        with open(path, "r+b") as f:
            f.truncate(size)
    return path


def _fasta(path, proteins=100, length=400):
    path.write_text("".join(f">sp|P{i}|X\n{'A' * length}\n" for i in range(proteins)))
    return path


def test_spectrum_count_from_header_or_size(tmp_path):
    assert spectrum_count(_mzml(tmp_path / "a.mzML", 1810)) == 1810
    no_count = tmp_path / "b.mzML"
    no_count.write_bytes(b"<mzML>" + b" " * (2 * 1024 * 1024))
    assert spectrum_count(no_count) == 240  # 2 MB at 120 spectra per MB


def test_fasta_stats(tmp_path):
    proteins, size = fasta_stats(_fasta(tmp_path / "db.fasta", proteins=7))
    assert proteins == 7 and size > 7 * 400


def test_search_factor():
    assert search_factor(REFERENCE) == 1.0
    assert search_factor({**REFERENCE, "num_enzyme_termini": "semi"}) == 8.0
    more_mods = {**REFERENCE, "variable_modifications": "Oxidation (M)\nAcetyl (Protein N-term)\nPhospho (S)\nPhospho (T)"}
    assert search_factor(more_mods) > 2.5
    open_search = {**REFERENCE, "precursor_mass_tolerance": "500", "precursor_error_units": "Da"}
    assert search_factor(open_search) > 100


def test_demo_sized_run_is_short_and_fits(tmp_path):
    files = [_mzml(tmp_path / f"{i}.mzML", 1810, 15 * 1024 * 1024) for i in range(2)]
    # the demo database: 3.7 MB, about 6,100 proteins
    est = estimate(files, _fasta(tmp_path / "db.fasta", 6100, 600), REFERENCE, max_threads=4, total_memory=16 * GB)
    assert est.parallel == 2 and est.spectra == 3620
    assert 0.5 < est.minutes < 5  # the demo run took about a minute
    assert 0.5 < est.peak_gb < 2  # Comet measured 550 MB per file, two at a time
    assert est.level == "ok" and not est.advice


def test_large_run_gets_memory_and_time_advice(tmp_path):
    files = [_mzml(tmp_path / f"{i}.mzML", 45000, 1024) for i in range(6)]
    fasta = _fasta(tmp_path / "db.fasta", 6100, 600)
    est = estimate(files, fasta, {**REFERENCE, "num_enzyme_termini": "semi"}, max_threads=4, total_memory=16 * GB)
    batch_mb = COMET_BASE_MB + COMET_MB_PER_SPECTRUM * 20000  # the default batch caps the spectra in memory
    assert abs(est.comet_gb_per_file - batch_mb / 1024) < 0.01
    assert est.level == "high"
    assert est.minutes > 120
    text = " ".join(est.advice)
    assert "spectrum_batch_size" in text and "fewer threads" in text and "two or three files" in text
    # a smaller batch brings it within memory
    smaller = estimate(files, fasta, {**REFERENCE, "spectrum_batch_size": "5000"}, max_threads=4, total_memory=16 * GB)
    assert smaller.peak_gb < 0.7 * 16


def test_format_minutes():
    assert format_minutes(0.4) == "about 1 minute"
    assert format_minutes(7) == "about 7 minutes"
    assert format_minutes(37) == "about 35 minutes"
    assert format_minutes(260) == "about 4.3 hours"
