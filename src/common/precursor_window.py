"""Shade the precursor XIC's m/z window on a PSM's spectrum plot.

The precursor XIC under the spectrum sums the MS1 signal within ``ppm`` of the
precursor m/z. LinePlot and MirrorPlot of openms-insight releases with
``window_identifier`` shade that window from the shared state; with an older
release nothing is passed and the plot is drawn as before.
"""
from __future__ import annotations

import inspect

PRECURSOR_WINDOW = "precursor_window"


def supports_window(component: type) -> bool:
    """Whether this openms-insight component takes ``window_identifier``."""
    try:
        return "window_identifier" in inspect.signature(component.__init__).parameters
    except (TypeError, ValueError):
        return False


def precursor_window(state_manager, mz: float | None, ppm: float, component: type) -> dict:
    """Arguments for ``component`` that shade ``mz`` ± ``ppm``; sets the window in the state.

    Gives ``{}`` (and leaves the state alone) without an m/z or when the
    installed openms-insight cannot shade a window.
    """
    if mz is None or not supports_window(component):
        return {}
    state_manager.set_selection(PRECURSOR_WINDOW, float(mz))
    return {"window_identifier": PRECURSOR_WINDOW, "window_tolerance": f"{ppm:g} ppm"}
