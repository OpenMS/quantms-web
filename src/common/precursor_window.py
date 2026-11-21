"""Shade a PSM's precursor isolation window on its spectrum plot, and mark its precursor m/z.

The isolation window is the m/z range the instrument selected for
fragmentation, an absolute range around the precursor: ``isolation window
target m/z`` minus its lower and plus its upper offset, as the mzML records
it per MS2 spectrum. LinePlot and MirrorPlot of openms-insight releases with
``window_identifier`` shade that range from the shared state; with an older
release nothing is passed and the plot is drawn as before.

Streamlit-free (pyOpenMS only).
"""
from __future__ import annotations

import inspect
import math

PRECURSOR_WINDOW = "precursor_window"


def _takes(component: type, option: str) -> bool:
    try:
        return option in inspect.signature(component.__init__).parameters
    except (TypeError, ValueError):
        return False


def supports_window(component: type) -> bool:
    """Whether this openms-insight component takes ``window_identifier``."""
    return _takes(component, "window_identifier")


def precursor_marker(mz, charge, component: type) -> dict:
    """Arguments for ``component`` that mark the precursor m/z with "precursor 2+".

    A dotted guide at the PSM's precursor m/z, labelled with its charge, drawn
    whether or not a peak is there and kept in the image download. Gives ``{}``
    without a usable m/z or when the installed openms-insight component has no
    ``x_markers``.
    """
    try:
        mz = float(mz)
    except (TypeError, ValueError):
        return {}
    if not math.isfinite(mz) or mz <= 0 or not _takes(component, "x_markers"):
        return {}
    try:
        label = f"precursor {int(charge)}+" if int(charge) > 0 else "precursor"
    except (TypeError, ValueError):
        label = "precursor"
    return {"x_markers": [{"x": mz, "label": label}]}


def isolation_windows(mzml_path) -> dict[int, tuple[float, float, float]]:
    """``{scan: (target, lower, upper)}``: each MS2 spectrum's isolation window in ``mzml_path``.

    ``target`` is the precursor m/z, ``lower`` and ``upper`` the offsets below
    and above it.

    Only spectrum headers are read (no peaks). Spectra without a precursor or
    without isolation offsets are left out.
    """
    import pyopenms as oms

    from src.common.identified_spectra import _scan

    options = oms.PeakFileOptions()
    options.setMSLevels([2])
    options.setFillData(False)
    exp = oms.MSExperiment()
    mzml = oms.MzMLFile()
    mzml.setOptions(options)
    mzml.load(str(mzml_path), exp)
    windows = {}
    for spec in exp:
        precursors = spec.getPrecursors()
        scan = _scan(spec.getNativeID())
        if not precursors or not scan:
            continue
        precursor = precursors[0]
        lower = precursor.getIsolationWindowLowerOffset()
        upper = precursor.getIsolationWindowUpperOffset()
        if lower > 0 or upper > 0:
            windows[scan] = (precursor.getMZ(), lower, upper)
    return windows


def precursor_window(state_manager, window: tuple[float, float, float] | None, component: type) -> dict:
    """Arguments for ``component`` that shade the isolation ``window``; sets it in the state.

    ``window`` is ``(target, lower, upper)``. A symmetric window is stored as
    its target with the offset as an m/z tolerance, so it reads
    ``652.8357 ± 0.75 m/z``; an asymmetric one as the range. Gives ``{}`` (and
    leaves the state alone) without a window or when the installed
    openms-insight cannot shade one.
    """
    if window is None or not supports_window(component):
        return {}
    target, lower, upper = (float(v) for v in window)
    if lower > 0 and math.isclose(lower, upper, abs_tol=1e-6):
        state_manager.set_selection(PRECURSOR_WINDOW, target)
        return {"window_identifier": PRECURSOR_WINDOW, "window_tolerance": f"{lower:g} m/z"}
    state_manager.set_selection(PRECURSOR_WINDOW, [target - lower, target + upper])
    return {"window_identifier": PRECURSOR_WINDOW}
