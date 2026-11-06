"""The instrument preset recommended from the selected runs (src/common/instrument.py)."""
from pathlib import Path

from src.common import instrument

DEMO = sorted(Path("example-data", "mzML").glob("*.mzML"))


def _mzml(tmp_path, name, model, analyzers, filters=()):
    components = "".join(f'<analyzer order="2"><cvParam cvRef="MS" accession="MS:1" name="{a}" /></analyzer>'
                         for a in analyzers)
    spectra = "".join(
        f'<spectrum><cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>'
        f'<cvParam cvRef="MS" accession="MS:1000512" name="filter string" value="{f}"/></spectrum>'
        for f in filters
    )
    text = (
        '<mzML><instrumentConfigurationList count="1"><instrumentConfiguration id="IC1">'
        f'<cvParam cvRef="MS" accession="MS:1002523" name="{model}" />'
        '<cvParam cvRef="MS" accession="MS:1000529" name="instrument serial number" value="1"/>'
        f'<componentList>{components}</componentList></instrumentConfiguration>'
        f'</instrumentConfigurationList><run>{spectra}</run></mzML>'
    )
    path = Path(tmp_path, name)
    path.write_text(text)
    return path


def test_demo_runs_are_orbitrap():
    rec = instrument.recommend([str(p) for p in DEMO])
    assert rec.preset == instrument.ORBITRAP
    assert rec.runs == len(DEMO)
    assert "Q Exactive HF" in rec.reason()


def test_ion_trap_ms2_wins_over_the_orbitrap_analyzer(tmp_path):
    path = _mzml(tmp_path, "velos.mzML", "LTQ Orbitrap Velos", ["orbitrap", "radial ejection linear ion trap"],
                 ["FTMS + p NSI Full ms [350.00-1800.00]", "ITMS + c NSI d Full ms2 500.25@cid35.00 [125.00-1500.00]"])
    run = instrument.read_mzml(path)
    assert run.preset == instrument.ION_TRAP
    assert run.model == "LTQ Orbitrap Velos"


def test_orbitrap_ms2_from_the_scan_filter(tmp_path):
    path = _mzml(tmp_path, "lumos.mzML", "Orbitrap Fusion Lumos", ["orbitrap", "radial ejection linear ion trap"],
                 ["FTMS + c NSI d Full ms2 500.25@hcd30.00 [120.00-1500.00]"])
    assert instrument.read_mzml(path).preset == instrument.ORBITRAP


def test_bruker_folder_is_timstof(tmp_path):
    folder = Path(tmp_path, "run.d")
    folder.mkdir()
    rec = instrument.recommend([str(folder)])
    assert rec.preset == instrument.TIMSTOF


def test_raw_counts_once_converted(tmp_path):
    raw = Path(tmp_path, "run.raw")
    raw.write_bytes(b"")
    converted = Path(tmp_path, "converted")
    rec = instrument.recommend([str(raw)], converted)
    assert rec.preset is None and rec.unreadable == 1
    converted.mkdir()
    _mzml(converted, "run.mzML", "Q Exactive", ["orbitrap"])
    rec = instrument.recommend([str(raw)], converted)
    assert rec.preset == instrument.ORBITRAP and rec.unreadable == 0


def test_mixed_instruments_recommend_nothing(tmp_path):
    folder = Path(tmp_path, "run.d")
    folder.mkdir()
    orbitrap = _mzml(tmp_path, "a.mzML", "Q Exactive", ["orbitrap"])
    rec = instrument.recommend([str(folder), str(orbitrap)])
    assert rec.preset is None
    assert rec.runs == 2
    assert "different instruments" in rec.reason()


def test_unknown_file_is_not_counted(tmp_path):
    path = Path(tmp_path, "x.mzML")
    path.write_text("<mzML></mzML>")
    rec = instrument.recommend([str(path), str(Path(tmp_path, "missing.mzML"))])
    assert rec.preset is None and rec.runs == 0 and rec.unreadable == 2


def test_every_preset_exists():
    import json
    presets = json.loads(Path("presets.json").read_text())["topp-workflow"]
    for name in instrument.PRESETS:
        assert {"ProSE", "ProSE-TMT"} <= set(presets[name])
