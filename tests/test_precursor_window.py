"""The precursor isolation window is shaded on the PSM spectrum when openms-insight can do it."""

import pytest

from src.common.precursor_window import PRECURSOR_WINDOW, isolation_windows, precursor_window, supports_window


class _State:
    def __init__(self):
        self.values = {}

    def set_selection(self, identifier, value):
        self.values[identifier] = value


class _WithWindow:
    def __init__(self, cache_id, window_identifier=None, window_tolerance=0.01, **kwargs):
        pass


class _WithoutWindow:
    def __init__(self, cache_id, **kwargs):
        pass


def test_symmetric_window_as_centre_and_mz_tolerance():
    state = _State()
    args = precursor_window(state, (825.74609375, 0.75, 0.75), _WithWindow)
    assert args == {"window_identifier": PRECURSOR_WINDOW, "window_tolerance": "0.75 m/z"}
    assert state.values == {PRECURSOR_WINDOW: pytest.approx(825.74609375)}


def test_asymmetric_window_as_range():
    state = _State()
    args = precursor_window(state, (612.3, 1.0, 0.5), _WithWindow)
    assert args == {"window_identifier": PRECURSOR_WINDOW}
    assert state.values == {PRECURSOR_WINDOW: pytest.approx([611.3, 612.8])}


def test_older_insight_or_no_window_gives_nothing():
    state = _State()
    assert precursor_window(state, (825.746, 0.75, 0.75), _WithoutWindow) == {}
    assert precursor_window(state, None, _WithWindow) == {}
    assert state.values == {}
    assert supports_window(_WithWindow) and not supports_window(_WithoutWindow)


def test_isolation_windows_from_mzml(tmp_path):
    oms = pytest.importorskip("pyopenms")
    exp = oms.MSExperiment()
    for scan, (mz, lower, upper) in {7: (825.74609375, 0.75, 0.75), 8: (500.25, 0.0, 0.0), 9: (612.3, 1.0, 0.5)}.items():
        spec = oms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
        precursor = oms.Precursor()
        precursor.setMZ(mz)
        precursor.setIsolationWindowLowerOffset(lower)
        precursor.setIsolationWindowUpperOffset(upper)
        spec.setPrecursors([precursor])
        exp.addSpectrum(spec)
    path = tmp_path / "run.mzML"
    oms.MzMLFile().store(str(path), exp)

    windows = isolation_windows(path)
    assert set(windows) == {7, 9}  # scan 8 records no isolation window
    assert windows[7] == pytest.approx((825.74609375, 0.75, 0.75))
    assert windows[9] == pytest.approx((612.3, 1.0, 0.5))
