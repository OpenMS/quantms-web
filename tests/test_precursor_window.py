"""The precursor XIC's m/z window is shaded on the PSM spectrum when openms-insight can do it."""

from src.common.precursor_window import PRECURSOR_WINDOW, precursor_window, supports_window


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


def test_window_set_and_passed():
    state = _State()
    args = precursor_window(state, 523.2711, 10.0, _WithWindow)
    assert args == {"window_identifier": PRECURSOR_WINDOW, "window_tolerance": "10 ppm"}
    assert state.values == {PRECURSOR_WINDOW: 523.2711}


def test_older_insight_or_no_mz_gives_nothing():
    state = _State()
    assert precursor_window(state, 523.2711, 10.0, _WithoutWindow) == {}
    assert precursor_window(state, None, 10.0, _WithWindow) == {}
    assert state.values == {}
    assert supports_window(_WithWindow) and not supports_window(_WithoutWindow)
