"""What the Viewer reads, computed at the end of a run.

Called by the workflow after its last tool, in the workflow process (no
Streamlit server; Streamlit's session state and caches work there in memory).
It fills the run's ``viewer-cache/`` (``src/common/disk_cache.py``) through the
same functions the Viewer calls, so the Viewer opens on finished results:
the per-file overview and QC tables, the per-sample counts, the protein,
peptide and PSM tables, the mass shifts, the abundance table and the group
comparison with the recommended settings.

Nothing here may stop a run: a failure is logged and the Viewer computes what
is missing when it is opened.
"""
from __future__ import annotations

from pathlib import Path


def _idxmls(directory: Path) -> list[Path]:
    return sorted(directory.glob("*.idXML")) if directory.is_dir() else []


def run(workspace: Path, run_id: str, log=print) -> None:
    """Compute and keep what the Viewer shows for ``run_id``."""
    import streamlit as st

    from src.common import overview, qc_plots, samples
    from src.common.guide import ensure_group_comparison
    from src.common.results_helpers import get_run_dir, get_run_params, search_results_dir
    from src.viewer.scope import experiment

    workspace = Path(workspace)
    run_dir = get_run_dir(workspace, run_id)
    results = run_dir / "results"
    params = get_run_params(workspace, run_id)
    mzml_files = [str(f) for f in params.get("mzML-files", []) or []]
    names = [Path(f).name for f in mzml_files]

    steps = [
        ("spectra", lambda: [overview.ms2_summary(f) for f in mzml_files]),
        ("identifications", lambda: [
            (overview.read_psms(p), samples.identification_sets(p)) for p in _idxmls(results / "psm_filter")
        ]),
        ("quality control", lambda: [
            qc_plots.psm_qc_table(p, names)
            for folder in (results / "percolator_results", search_results_dir(results), results / "psm_filter")
            for p in _idxmls(folder)
        ]),
    ]

    def comparison():
        previous = st.session_state.get("workspace")
        st.session_state["workspace"] = workspace
        try:
            # The run on its own: the comparison its Viewer opens on
            with experiment(run_id, 0, single=True):
                ensure_group_comparison(params, show_errors=False)
        finally:
            if previous is None:
                st.session_state.pop("workspace", None)
            else:
                st.session_state["workspace"] = previous

    def in_viewer_scope(compute):
        """``compute()`` with the workspace and the run picked as the Viewer has them."""
        def step():
            previous = st.session_state.get("workspace")
            st.session_state["workspace"] = workspace
            try:
                with experiment(run_id, 0, single=True):
                    return compute()
            finally:
                if previous is None:
                    st.session_state.pop("workspace", None)
                else:
                    st.session_state["workspace"] = previous
        return step

    def protein_tables():
        from types import SimpleNamespace

        from src.viewer import protein_data

        protein_data.load(SimpleNamespace(workspace=workspace))

    def mass_shifts():
        from src.viewer.panels import mass_shifts as panel

        idxmls = _idxmls(results / "psm_filter")
        if idxmls:
            panel.summary(tuple(str(p) for p in idxmls), tuple(names), tuple(p.stat().st_mtime for p in idxmls))

    steps.append(("protein tables", in_viewer_scope(protein_tables)))
    steps.append(("mass shifts", mass_shifts))
    steps.append(("group comparison", comparison))
    for label, step in steps:
        try:
            step()
        except Exception as err:  # pylint: disable=broad-except
            log(f"Viewer: {label} not prepared ({err}); it is computed when the Viewer opens.")
