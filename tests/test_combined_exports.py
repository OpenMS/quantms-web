"""Tests for src/common/combined_exports.py: the merged idXML and design.tsv of a run."""
import json
import shutil
from pathlib import Path

import pandas as pd
import pyopenms as poms
import pytest

from src.common.combined_exports import (
    build_missing,
    combined_paths,
    isobaric_design_tsv,
    merge_command,
    run_design_tsv,
    write_design,
)
from src.workflow.ExperimentalDesign import read_design

FILES = [
    "../ws/input-files/mzML-files/a_f1.mzML",
    "../ws/input-files/mzML-files/a_f2.mzML",
    "../ws/input-files/mzML-files/b_f1.mzML",
    "../ws/input-files/mzML-files/b_f2.mzML",
]


def _lfq_params():
    return {
        "analysis-mode": "LFQ",
        "mzML-files": FILES,
        "mzML-group-a_f1.mzML": "ctrl",
        "mzML-group-a_f2.mzML": "ctrl",
        "mzML-group-b_f1.mzML": "treated",
        "mzML-group-b_f2.mzML": "treated",
    }


def _fractions():
    return pd.DataFrame({
        "file": ["a_f1.mzML", "a_f2.mzML", "b_f1.mzML", "b_f2.mzML"],
        "condition": ["old", "old", "old", "old"],
        "replicate": [1, 1, 1, 1],
        "fraction": [1, 2, 1, 2],
    })


def test_lfq_design_uses_run_conditions_and_bare_names(tmp_path):
    path = write_design(tmp_path, _lfq_params(), _fractions())
    assert path == combined_paths(tmp_path)["design"]
    text = path.read_text()
    assert "../ws" not in text  # portable: file names only
    design = read_design(path)
    assert list(design["condition"]) == ["ctrl", "ctrl", "treated", "treated"]
    assert list(design["fraction"]) == [1, 2, 1, 2]
    ed = poms.ExperimentalDesignFile.load(str(path), False)
    assert ed.getNumberOfSamples() == 2
    assert ed.getNumberOfFractions() == 2


def test_lfq_design_without_stored_table_is_one_sample_per_file(tmp_path):
    path = write_design(tmp_path, _lfq_params())
    ed = poms.ExperimentalDesignFile.load(str(path), False)
    assert ed.getNumberOfSamples() == 4
    assert ed.getNumberOfFractions() == 1


def test_isobaric_design_one_sample_per_file_and_channel(tmp_path):
    params = {"TMT-group-sample1": "ctrl", "TMT-group-sample2": "skip", "TMT-group-sample3": "treated"}
    text = isobaric_design_tsv(FILES[:2], params, "tmt6plex")
    path = tmp_path / "design.tsv"
    path.write_text(text)
    ed = poms.ExperimentalDesignFile.load(str(path), False)
    assert ed.getNumberOfSamples() == 12
    assert ed.getNumberOfLabels() == 6
    samples = text.split("\n\n")[1]
    rows = [line.split("\t") for line in samples.strip().split("\n")[1:]]
    assert [r[1] for r in rows[:4]] == ["ctrl", "NA", "treated", "NA"]


def test_run_design_tsv_picks_isobaric_for_tmt_runs():
    params = {"analysis-mode": "TMT", "mzML-files": FILES[:1], "IsobaricAnalyzer-TMT": {"type": "tmt10plex"}}
    file_section = run_design_tsv(params).split("\n\n")[0]
    assert len(file_section.strip().split("\n")) == 1 + 10


def test_merge_command_keeps_file_origin():
    cmd = merge_command(["x.idXML", "y.idXML"], "out.idXML")
    assert cmd[0] == "IDMerger"
    assert cmd[cmd.index("-in") + 1:cmd.index("-out")] == ["x.idXML", "y.idXML"]
    assert cmd[cmd.index("-annotate_file_origin") + 1] == "true"


def test_build_missing_creates_both_for_an_older_run(tmp_path):
    results = tmp_path / "results"
    (results / "psm_filter").mkdir(parents=True)
    for name in ("b_filter.idXML", "a_filter.idXML"):
        (results / "psm_filter" / name).write_text("<IdXML/>")
    calls = []

    def fake_run(cmd):
        calls.append(cmd)
        Path(cmd[cmd.index("-out") + 1]).write_text("merged")
        return True

    done = build_missing(results, _lfq_params(), tmp_path, run=fake_run)
    assert done == {"idxml": True, "design": True}
    # files are merged in a stable (sorted) order
    assert [Path(p).name for p in calls[0][2:4]] == ["a_filter.idXML", "b_filter.idXML"]
    assert combined_paths(results)["design"].is_file()
    # nothing left to do the second time
    assert build_missing(results, _lfq_params(), tmp_path, run=fake_run) == {}
    assert len(calls) == 1


def test_build_missing_prefers_the_design_proteomicslfq_got(tmp_path):
    results = tmp_path / "results"
    (results / "quant_results").mkdir(parents=True)
    write_design(tmp_path / "x", _lfq_params(), _fractions())
    shutil.copy(tmp_path / "x" / "combined" / "design.tsv", results / "quant_results" / "experimental-design.tsv")
    build_missing(results, _lfq_params(), tmp_path, run=lambda c: False)
    assert list(read_design(combined_paths(results)["design"])["fraction"]) == [1, 2, 1, 2]


@pytest.mark.skipif(shutil.which("IDMerger") is None, reason="needs the IDMerger TOPP tool")
def test_idmerger_merges_two_runs(tmp_path):
    for i, name in enumerate(("a", "b")):
        prot = poms.ProteinIdentification()
        prot.setIdentifier(f"run{i}")
        pep = poms.PeptideIdentification()
        pep.setIdentifier(f"run{i}")
        hit = poms.PeptideHit()
        hit.setSequence(poms.AASequence.fromString("PEPTIDEK"))
        hit.setCharge(2)
        pep.setHits([hit])
        peps = poms.PeptideIdentificationList() if hasattr(poms, "PeptideIdentificationList") else []
        peps.push_back(pep) if hasattr(peps, "push_back") else peps.append(pep)
        poms.IdXMLFile().store(str(tmp_path / f"{name}.idXML"), [prot], peps)
    out = tmp_path / "merged.idXML"
    from src.common.combined_exports import _run
    assert _run(merge_command([tmp_path / "a.idXML", tmp_path / "b.idXML"], out))
    prots = []
    peps = poms.PeptideIdentificationList() if hasattr(poms, "PeptideIdentificationList") else []
    poms.IdXMLFile().load(str(out), prots, peps)
    assert len(peps) == 2
    assert json.dumps(sorted(str(p.getMetaValue("file_origin")) for p in peps)).count(".idXML") == 2
