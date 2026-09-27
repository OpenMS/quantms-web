"""Tests for the PSM / fragment ion exports and the streamed result tables."""
from pathlib import Path

import pandas as pd
import polars as pl
import pytest

openms_insight = pytest.importorskip("openms_insight")
from openms_insight import SequenceView  # noqa: E402

if not hasattr(SequenceView, "export_fragment_ions"):
    pytest.skip("needs openms-insight with SequenceView.export_fragment_ions", allow_module_level=True)

from src.common.results_helpers import write_psm_exports  # noqa: E402

PROTON = 1.007276


def _y_ion_mz(sequence: str, n: int, charge: int = 1) -> float:
    import pyopenms as poms

    neutral = poms.AASequence.fromString(sequence).getSuffix(n).getMonoWeight()
    return (neutral + charge * PROTON) / charge


def test_write_psm_exports(tmp_path):
    id_df = pl.DataFrame(
        {
            "id_idx": [0, 1],
            "scan_id": [5, 6],
            "file_index": [0, 1],
            "filename": ["A.mzML", "B.mzML"],
            "sequence": ["PEPTIDEK", "ELVISK"],
            "charge": [2, 2],
            "mz": [465.2, 351.7],
            "rt": [100.0, 200.0],
            "score": [0.001, 0.002],
            "protein_accession": ["P1", "P2"],
        }
    )
    peaks = pl.DataFrame(
        {
            "peak_id": [0, 1, 2],
            "file_index": [0, 1, 0],
            "scan_id": [5, 6, 6],
            "mass": [_y_ion_mz("PEPTIDEK", 3), _y_ion_mz("ELVISK", 2), 999.0],
            "intensity": [10.0, 20.0, 30.0],
        }
    )
    seq_view = SequenceView(
        cache_id="seqview_test",
        sequence_data=id_df.lazy()
        .select(["id_idx", "sequence", "charge", "file_index", "scan_id"])
        .rename({"id_idx": "sequence_id", "charge": "precursor_charge"}),
        peaks_data=peaks.lazy(),
        filters={"identification": "sequence_id", "file": "file_index", "spectrum": "scan_id"},
        cache_path=str(tmp_path / "cache"),
        deconvolved=False,
        annotation_config={"ion_types": ["b", "y"], "tolerance": 0.02, "tolerance_ppm": False},
    )

    write_psm_exports(id_df, seq_view, tmp_path / "exports", "run_filter")

    psms = pl.read_csv(tmp_path / "exports" / "run_filter_psms.tsv", separator="\t")
    assert psms.height == 2
    assert set(id_df.columns) <= set(psms.columns)

    frags = pl.read_csv(tmp_path / "exports" / "run_filter_fragment_ions.tsv", separator="\t")
    matched = {(r["sequence"], r["ion"], r["filename"]) for r in frags.iter_rows(named=True) if r["charge"] == 1}
    assert ("PEPTIDEK", "y3", "A.mzML") in matched
    assert ("ELVISK", "y2", "B.mzML") in matched
    # PSM context comes along with every ion
    assert {"rt", "precursor_mz", "score", "protein_accession"} <= set(frags.columns)


def test_insight_table_cache_follows_source_file(tmp_path, monkeypatch):
    from src.common import insight_tables

    rendered = []
    monkeypatch.setattr(insight_tables.Table, "__call__", lambda self, **kw: rendered.append(self._cache_id))

    source = tmp_path / "quant.csv"
    source.write_text("x")
    df = pd.DataFrame({"Protein.Name": ["A", "B"], "sample1": [1.0, 2.0], "note": [None, "x"]})

    insight_tables.show_insight_table(df, "t", tmp_path / "cache", source_files=[source])
    first = rendered[-1]
    table = insight_tables.Table(cache_id=first, cache_path=str(tmp_path / "cache"))
    fields = [c["field"] for c in table._column_definitions]
    assert fields == ["Protein_Name", "sample1", "note"]
    assert table._column_definitions[0]["title"] == "Protein.Name"

    # A rerun rewrites the source: a new cache replaces the old one
    import os
    os.utime(source, ns=(1, 1))
    insight_tables.show_insight_table(df, "t", tmp_path / "cache", source_files=[source])
    assert rendered[-1] != first
    assert not (tmp_path / "cache" / first).exists()
