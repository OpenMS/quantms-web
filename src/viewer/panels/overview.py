"""results_overview as a Viewer component."""

from pathlib import Path
import streamlit as st
from src.common.empty_result import fdr_summary
from src.common.insight_tables import show_insight_table
from src.common.id_diagnosis import EXPECTED_ID_RATE, diagnose_identifications, search_kind
from src.common.overview import DIAGNOSIS_COLUMNS, run_overview
from src.common.results_helpers import get_abundance_data, get_run_dir, get_run_params, search_tolerances
from src.viewer.scope import needs_results


def render(ctx) -> None:
    run_dir = get_run_dir(ctx.workspace)
    filter_dir = run_dir / "results" / "psm_filter"
    cache_dir = run_dir / "results" / "insight_cache"
    workflow_params = get_run_params(ctx.workspace)
    mzml_files = [str(f) for f in workflow_params.get("mzML-files", [])]
    if not filter_dir.is_dir() or not any(filter_dir.glob("*.idXML")):
        return needs_results()
    log_file = run_dir / "logs" / "commands-and-run-times.log"
    sources = sorted(filter_dir.glob("*.idXML")) + [log_file]

    def _mtimes(paths) -> tuple:
        return tuple(Path(p).stat().st_mtime if Path(p).exists() else 0.0 for p in paths)

    tolerances = search_tolerances(workflow_params)
    precursor_ppm = tolerances["precursor"] if tolerances["precursor_unit"] == "ppm" else None

    @st.cache_data(show_spinner="Summarizing the run...")
    def _overview(workflow_dir: str, files: tuple, precursor_ppm, mtimes: tuple):
        """run_overview, cached on the inputs' and results' modification times."""
        return run_overview(workflow_dir, list(files), precursor_ppm=precursor_ppm)

    per_file, figures = _overview(
        str(run_dir), tuple(mzml_files), precursor_ppm, _mtimes(mzml_files + [str(p) for p in sources])
    )
    abundance = get_abundance_data(ctx.workspace)
    proteins = samples = conditions = None
    if abundance is not None:
        pivot_df, expr_df, group_map = abundance
        proteins = len(pivot_df)
        samples = len(expr_df.columns)
        conditions = len({str(g) for g in group_map.values() if str(g).strip()}) or None
    fdr = float(workflow_params.get("psm-fdr-percent", 1.0))

    def _show(value, unit: str = "") -> str:
        return "–" if value is None else f"{value:,}{unit}"

    st.subheader("Key figures")
    st.caption(f"Filters of this run: {fdr_summary(workflow_params)}.")
    row1 = st.columns(4)
    row1[0].metric("Quantified proteins", _show(proteins), help="Proteins in the quantification result.")
    row1[1].metric(
        "Peptides",
        _show(figures["Peptides"]),
        help=f"Distinct target peptide sequences at {fdr:g} % PSM FDR, over all files.",
    )
    row1[2].metric(
        "Target PSMs", _show(figures["Target PSMs"]), help=f"Target PSMs at {fdr:g} % PSM FDR, over all files."
    )
    row1[3].metric(
        "PSMs per MS2 spectrum",
        _show(figures["PSMs per MS2 (%)"], " %"),
        help="Target PSMs divided by the MS2 spectra of all files.",
    )
    row2 = st.columns(4)
    row2[0].metric(
        "Median precursor error",
        _show(figures["Median precursor error (ppm)"], " ppm"),
        help="Median m/z error of the target PSMs, after correcting picks of a 13C isotope peak.",
    )
    row2[1].metric("Files", _show(figures["Files"]))
    row2[2].metric("Samples", _show(samples), help="Sample columns of the quantification result.")
    row2[3].metric("Conditions", _show(conditions), help="Distinct groups set on the Describe experiment page.")
    _diagnosis(per_file, workflow_params, cache_dir, sources)
    st.subheader("Files")
    st.caption(
        f"PSM counts are at {fdr:g} % PSM FDR. ID rate is target PSMs per MS2 spectrum. Mass error is "
        "the median precursor m/z error of the file's target PSMs. Contaminant PSMs is the share of "
        "target PSMs on contaminants other than trypsin. Percolator is whether it rescored "
        "the file. Search is the database search's (ProSE) run time in seconds in the last run."
    )
    # One row per file; the fractional columns go as text, since the table stores
    # floats as float32 and would show 36.8 as 36.79999923706055
    per_file_view = per_file.drop(columns=[c for c in DIAGNOSIS_COLUMNS if c != "Contaminant PSMs (%)"])
    for col in ("ID rate (%)", "Mass error (ppm)", "Contaminant PSMs (%)"):
        per_file_view[col] = per_file_view[col].map(lambda v: "" if v is None or v != v else f"{v:g}")
    show_insight_table(
        per_file_view,
        name="overview_files",
        cache_dir=cache_dir,
        source_files=sources,
        title="run_overview",
        height=300,
        # Share the page width between the columns rather than folding the last ones away
        layout="fitColumns",
    )


def _diagnosis(per_file, workflow_params: dict, cache_dir, sources) -> None:
    """The files whose ID rate falls short and their likely causes; nothing when none does."""
    flagged = diagnose_identifications(per_file, workflow_params)
    if flagged.empty:
        return
    kind = search_kind(workflow_params)
    st.subheader("Why so few identifications?")
    st.caption(
        f"{len(flagged)} of {len(per_file)} files fall short: below the {EXPECTED_ID_RATE[kind]:g} % "
        f"a {kind} search usually reaches, or below a third of the run's median. The likely causes "
        "come from the files' spectra and PSMs, the most likely first."
    )
    view = flagged.copy()
    view["ID rate (%)"] = view["ID rate (%)"].map(lambda v: "" if v is None or v != v else f"{v:g}")
    view["Expected (%)"] = view["Expected (%)"].map(lambda v: f"{v:g}")
    show_insight_table(
        view,
        name="overview_diagnosis",
        cache_dir=cache_dir,
        source_files=sources,
        extra_key=kind,
        title="id_diagnosis",
        # The causes are sentences: wrap them in a wide column rather than cut them off
        height=min(600, 90 + 95 * len(view)),
        layout="fitColumns",
        column_overrides={
            "File": {"widthGrow": 2, "hozAlign": "left"},
            "ID rate (%)": {"widthGrow": 1},
            "Expected (%)": {"widthGrow": 1},
            "Likely causes": {"formatter": "textarea", "widthGrow": 7, "hozAlign": "left"},
        },
    )
