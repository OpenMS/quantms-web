"""Helper functions for results pages."""
import json
import re
import pandas as pd
import polars as pl
import numpy as np
import streamlit as st
from pathlib import Path
from pyopenms import IdXMLFile, MSExperiment, MzMLFile
try:  # pyOpenMS >= 3.5.0
    from pyopenms import PeptideIdentificationList
except ImportError:  # pyOpenMS <= 3.4.x
    PeptideIdentificationList = None
from src.workflow.ParameterManager import ParameterManager
from src.workflow.RunManager import RunManager
from src.workflow.run_widgets import get_selected_run, select_run
from src.common.contaminants import drop_contaminants, is_contaminant
from src.common.peptide_format import to_proforma
from src.common.isobaric import channel_column_names, label_type

UNASSIGNED = "Unassigned"
# Columns the group comparison adds to a protein table
STAT_COLUMNS = ("log2FC", "p-value", "p-adj", "stat")

def get_workflow_dir(workspace):
    """Get the workflow directory path."""
    return Path(workspace, "topp-workflow")


# Session-state keys holding results computed from the shown run. They are
# dropped when another run is picked, so no page mixes two runs.
RUN_DERIVED_KEYS = [
    "filtered_df", "imputed_df", "normalized_df", "statistics_df",
    "normalized_transform", "comparison_settings", "comparison_group_map",
    "go_enrichment", "go_results", "go_ready", "pivot_df",
]
# Group assignments are edited on Describe experiment after a run, so they are
# read from the current settings; every other parameter as the run was started.
GROUP_PREFIXES = ("mzML-group-", "TMT-group-")


def get_runs(workspace) -> RunManager:
    """The analysis runs of the workspace (a results/ from before runs became one)."""
    runs = RunManager(get_workflow_dir(workspace))
    runs.adopt_legacy_results()
    return runs


def get_run_id(workspace) -> str | None:
    """The run the results pages show (the newest finished one unless picked)."""
    return get_selected_run(get_runs(workspace))


def get_results_dir(workspace, run_id: str | None = None) -> Path:
    """Results folder of a run (default: the one shown). Missing if nothing ran yet."""
    runs = get_runs(workspace)
    run_id = run_id or get_selected_run(runs)
    if run_id is None:
        return get_workflow_dir(workspace) / "results"
    return runs.results_dir(run_id)


def get_run_dir(workspace, run_id: str | None = None) -> Path:
    """Folder of a run (default: the one shown): results/, logs/ and params.json,
    laid out like the workflow folder."""
    return get_results_dir(workspace, run_id).parent


def get_run_params(workspace, run_id: str | None = None) -> dict:
    """Parameters of a run (default: the one shown), with the current sample groups."""
    workflow_dir = get_workflow_dir(workspace)
    current = ParameterManager(workflow_dir, "TOPP Workflow").get_parameters_from_json()
    runs = get_runs(workspace)
    run_id = run_id or get_selected_run(runs)
    params = runs.params(run_id) if run_id else {}
    if not params:
        return current
    return params | {k: v for k, v in current.items() if k.startswith(GROUP_PREFIXES)}


def run_picker(workspace=None) -> str | None:
    """Sidebar picker for the run the results pages show. Returns its id.

    Results derived from another run (the group comparison, GO enrichment) are
    dropped when the pick changes, including when a new run has finished.
    """
    workspace = workspace or st.session_state.get("workspace")
    if workspace is None:
        return None
    runs = get_runs(workspace)
    with st.sidebar:
        run_id = select_run(
            runs, label="Showing run",
            help="Every analysis you start is kept as its own run. Pick the one the "
            "results pages show; rename, download or delete runs on Runs and downloads.",
        )
    if st.session_state.get("_results_run") != run_id:
        for key in RUN_DERIVED_KEYS:
            st.session_state.pop(key, None)
        st.session_state["_results_run"] = run_id
    return run_id


def load_idxml(idxml_file):
    """Load an idXML file, returning ``(protein_ids, peptide_ids)`` as plain lists.

    pyOpenMS 3.5.0 changed the third parameter of ``IdXMLFile.load()`` from a
    ``libcpp_vector[PeptideIdentification]``, which accepted an ordinary Python
    list, to a dedicated ``PeptideIdentificationList`` container. A list no longer
    matches any overload there and pyOpenMS raises
    ``Exception: can not handle type of (<path>, [], [])``. Feature-detect the
    container so the app runs on 3.5.0 as well as earlier releases, and hand
    callers back an ordinary list either way.
    """
    proteins = []
    peptides = PeptideIdentificationList() if PeptideIdentificationList else []
    IdXMLFile().load(str(idxml_file), proteins, peptides)
    return proteins, list(peptides)


def idxml_to_df(idxml_file):
    """Parse idXML file and return DataFrame with peptide hits."""
    proteins, peptides = load_idxml(idxml_file)

    records = []
    for pep in peptides:
        rt = pep.getRT()
        mz = pep.getMZ()
        for h in pep.getHits():
            protein_refs = [ev.getProteinAccession() for ev in h.getPeptideEvidences()]
            records.append({
                "RT": rt,
                "m/z": mz,
                "Sequence": h.getSequence().toString(),
                "Charge": h.getCharge(),
                "Score": h.getScore(),
                "Proteins": ",".join(protein_refs) if protein_refs else None,
            })

    df = pd.DataFrame(records)
    if not df.empty:
        df["Charge"] = df["Charge"].astype(str)
        df["Charge_num"] = df["Charge"].astype(int)
    return df


def create_psm_scatter_plot(df_plot):
    """Create a scatter plot for PSM visualization."""
    import plotly.express as px

    fig = px.scatter(
        df_plot,
        x="RT",
        y="m/z",
        color="Score",
        custom_data=["index", "Sequence", "Proteins"],
        color_continuous_scale=["#a6cee3", "#1f78b4", "#08519c", "#08306b"],
    )
    fig.update_traces(
        marker=dict(size=6, opacity=0.8),
        hovertemplate='<b>Index: %{customdata[0]}</b><br>'
                    + 'RT: %{x:.2f}<br>'
                    + 'm/z: %{y:.4f}<br>'
                    + 'Score: %{marker.color:.3f}<br>'
                    + 'Sequence: %{customdata[1]}<br>'
                    + 'Proteins: %{customdata[2]}<br>'
                    + '<extra></extra>'
    )
    fig.update_layout(
        coloraxis_colorbar=dict(title="Score"),
        hovermode="closest"
    )
    return fig


def extract_scan_from_ref(spec_ref: str) -> int:
    """Extract scan number from spectrum reference string.

    Format: "controllerType=0 controllerNumber=1 scan=1234"
    """
    match = re.search(r'scan=(\d+)', spec_ref)
    return int(match.group(1)) if match else 0


def extract_scan_number(native_id: str) -> int:
    """Extract scan number from native ID."""
    match = re.search(r'scan=(\d+)', native_id)
    return int(match.group(1)) if match else 0


def extract_filename_from_idxml(idxml_path: Path, mzml_names=None) -> str:
    """Derive the source mzML filename from an idXML filename.

    The workflow names idXML files ``<mzML stem>_<step>...``, e.g.
    ``02COVID_comet_perc_filter.idXML``. When the input mzML names are known the
    longest one whose stem prefixes the idXML stem wins; otherwise the known
    step suffixes are stripped from the end.
    """
    stem = idxml_path.stem
    if mzml_names:
        candidates = [
            name for name in mzml_names
            if stem == Path(name).stem or stem.startswith(f"{Path(name).stem}_")
        ]
        if candidates:
            return max(candidates, key=len)
    while True:
        for suffix in ['_filter', '_perc', '_per', '_comet']:
            if stem.endswith(suffix):
                stem = stem[:-len(suffix)]
                break
        else:
            return f"{stem}.mzML"


ID_DF_SCHEMA = {
    "id_idx": pl.Int64,
    "scan_id": pl.Int64,
    "file_index": pl.Int64,
    "filename": pl.Utf8,
    "sequence": pl.Utf8,
    "charge": pl.Int64,
    "mz": pl.Float64,
    "rt": pl.Float64,
    "score": pl.Float64,
    "protein_accession": pl.Utf8,
    "contaminant": pl.Boolean,
}


def parse_idxml(idxml_path: Path, filename_to_index: dict | None = None) -> tuple[pl.DataFrame, list[str]]:
    """Parse idXML and return DataFrame for openms_insight.

    Args:
        idxml_path: idXML file of one run.
        filename_to_index: Maps every input mzML filename to its file index in
            the spectra cache. The PSMs get the index of their own run, so they
            link to that run's spectra.

    Returns:
        Tuple of (id_df, spectra_data list of source filenames)
    """
    proteins, peptides = load_idxml(idxml_path)

    # Derive mzML filename from idXML filename (e.g., 02COVID_filter.idXML -> 02COVID.mzML)
    source = extract_filename_from_idxml(idxml_path, list(filename_to_index or []))
    spectra_data = [source]

    run_index = (filename_to_index or {}).get(source, 0)

    records = []
    for pep in peptides:
        # Get spectrum reference from meta value (key may be bytes or string)
        spec_ref = ""
        if pep.metaValueExists("spectrum_reference"):
            spec_ref = pep.getMetaValue("spectrum_reference")
            if isinstance(spec_ref, bytes):
                spec_ref = spec_ref.decode()
        scan_id = extract_scan_from_ref(spec_ref)

        # One idXML holds one run, so every PSM belongs to that run's file index
        file_index = run_index
        filename = spectra_data[0]

        for h in pep.getHits():
            accessions = ";".join([ev.getProteinAccession() for ev in h.getPeptideEvidences()])
            records.append({
                "id_idx": len(records),
                "scan_id": scan_id,
                "file_index": file_index,
                "filename": Path(filename).name if filename else "",
                "sequence": h.getSequence().toString(),
                "charge": h.getCharge(),
                "mz": pep.getMZ(),
                "rt": pep.getRT(),
                "score": h.getScore(),
                "protein_accession": accessions,
                "contaminant": is_contaminant(accessions),
            })

    # Keep the columns when nothing was identified, so callers selecting
    # them get an empty table rather than a ColumnNotFoundError.
    return pl.DataFrame(records, schema=ID_DF_SCHEMA), spectra_data


def build_spectra_cache(mzml_dir: Path, filename_to_index: dict, mzml_files=None) -> tuple[pl.DataFrame, dict]:
    """Extract MS2 spectra from mzML files and return DataFrame.

    Args:
        mzml_dir: Directory containing mzML files
        filename_to_index: Dict mapping filename to file_index
        mzml_files: mzML paths to read. Defaults to every mzML in ``mzml_dir``;
            pass the workflow inputs so referenced (not copied) files are found.

    Returns:
        Tuple of (spectra_df, updated filename_to_index)
    """
    records = []
    peak_id = 0

    paths = [Path(f) for f in mzml_files] if mzml_files else sorted(mzml_dir.glob("*.mzML"))
    for mzml_path in paths:
        # Get or create file index
        if mzml_path.name not in filename_to_index:
            filename_to_index[mzml_path.name] = len(filename_to_index)
        file_index = filename_to_index[mzml_path.name]

        exp = MSExperiment()
        MzMLFile().load(str(mzml_path), exp)

        for spec in exp:
            if spec.getMSLevel() != 2:
                continue
            scan_id = extract_scan_number(spec.getNativeID())
            mz_array, int_array = spec.get_peaks()

            for mz, intensity in zip(mz_array, int_array):
                records.append({
                    "peak_id": peak_id,
                    "file_index": file_index,
                    "scan_id": scan_id,
                    "mass": float(mz),      # Use "mass" not "mz"
                    "intensity": float(intensity),
                })
                peak_id += 1

    return pl.DataFrame(records), filename_to_index


def write_psm_exports(id_df: pl.DataFrame, seq_view, export_dir: Path, prefix: str) -> None:
    """Write one run's PSM table and per-scan fragment ion table as TSV.

    ``<prefix>_psms.tsv`` holds every PSM with all its columns.
    ``<prefix>_fragment_ions.tsv`` holds, for every PSM, each fragment ion the
    SequenceView matched in its spectrum (same ion types, losses and tolerance
    as the view), together with the PSM's run, RT, m/z, score and proteins.
    """
    export_dir.mkdir(parents=True, exist_ok=True)
    # Exports carry the ProForma form next to the OpenMS string, for other tools
    if "sequence" in id_df.columns:
        id_df = id_df.with_columns(
            pl.col("sequence").map_elements(to_proforma, return_dtype=pl.Utf8).alias("proforma")
        )
    id_df.write_csv(export_dir / f"{prefix}_psms.tsv", separator="\t")

    fragments = seq_view.export_fragment_ions()
    if fragments.height and "sequence_id" in fragments.columns and "id_idx" in id_df.columns:
        psm_info = id_df.select(
            pl.col("id_idx").cast(fragments.schema["sequence_id"]).alias("sequence_id"),
            "filename",
            "rt",
            "proforma",
            pl.col("mz").alias("precursor_mz"),
            "score",
            "protein_accession",
        )
        fragments = fragments.join(psm_info, on="sequence_id", how="left")
    fragments.write_csv(export_dir / f"{prefix}_fragment_ions.tsv", separator="\t")


@st.cache_data
def load_abundance_data(
    workspace_path: str, csv_mtime: float, params_mtime: float = 0.0, run_id: str | None = None
) -> tuple | None:
    """Load CSV and build abundance matrices for downstream preprocessing.

    Args:
        workspace_path: Path to the workspace directory
        csv_mtime: Modification time of CSV file (used as cache key)
        params_mtime: Modification time of params.json (used as cache key so
            changing group assignments in Configure invalidates the cache
            even when the CSV itself hasn't changed)
        run_id: The run to load (default: the newest finished one)

    Returns:
        Tuple of (pivot_df, expr_df, group_map) or None if data unavailable
    """
    quant_dir = get_results_dir(workspace_path, run_id) / "quant_results"
    workflow_params = get_run_params(workspace_path, run_id)
    analysis_mode = workflow_params.get("analysis-mode", "LFQ")

    if analysis_mode == "LFQ":
        if not quant_dir.exists():
            return None

        csv_files = sorted(quant_dir.glob("*.csv"))
        if not csv_files:
            return None

        csv_file = csv_files[0]

        try:
            df = pd.read_csv(csv_file)
        except Exception:
            return None

        if df.empty:
            return None

        # Get optional group mapping from parameters.
        # Group information is not required at this stage; statistical testing
        # happens in the Statistical page.
        params = workflow_params
        group_map = {
            key[11:]: value  # Remove "mzML-group-" prefix
            for key, value in params.items()
            if key.startswith("mzML-group-") and value
        }

        df["Sample"] = df["Reference"].str.replace(".mzML", "", regex=False)

        # Build sample display order.
        if group_map:
            sample_group_df = df[["Sample", "Reference"]].drop_duplicates()
            sample_group_df["Group"] = sample_group_df["Reference"].map(group_map)
            grouped_samples = []
            for grp in sorted(sample_group_df["Group"].dropna().unique()):
                grouped_samples.extend(
                    sample_group_df[sample_group_df["Group"] == grp]["Sample"].tolist()
                )
            remaining_samples = [
                s for s in sorted(df["Sample"].unique()) if s not in grouped_samples
            ]
            all_samples = grouped_samples + remaining_samples
        else:
            all_samples = sorted(df["Sample"].unique())

        # Build pivot table
        pivot_list = []
        for protein, group_df in df.groupby("ProteinName"):
            peptides = ";".join(group_df["PeptideSequence"].unique())
            intensity_dict = group_df.groupby("Sample")["Intensity"].sum().to_dict()
            intensity_dict_complete = {
                sample: intensity_dict.get(sample, 0)
                for sample in all_samples
            }
            row = {
                "ProteinName": protein,
                **intensity_dict_complete,
                "PeptideSequence": peptides,
            }
            pivot_list.append(row)

        pivot_df = pd.DataFrame(pivot_list)
        pivot_df = pivot_df[["ProteinName"] + all_samples + ["PeptideSequence"]]

        # Build expression matrix (log2-transformed)
        expr_df = pivot_df.set_index("ProteinName")[all_samples]
        expr_df = expr_df.replace(0, np.nan)
        expr_df = np.log2(expr_df + 1)
        expr_df = expr_df.dropna()

        return pivot_df, expr_df, group_map

    else:
        if not quant_dir.exists():
            return None

        csv_files = sorted(quant_dir.glob("*.csv"))
        if not csv_files:
            return None

        csv_file = csv_files[0]

        try:
            df = pd.read_csv(csv_file, sep="\t", comment="#", engine="python")
        except Exception:
            return None

        if df.empty:
            return None

        # ratio column removal
        df = df.loc[:, ~df.columns.str.contains('ratio', case=False)]
        # "abundance_sample3[unknown]" -> "sample3[127C]"
        df = df.rename(columns=channel_column_names(df.columns, label_type(workflow_params)))

        # exclude_indices = st.session_state.get("tmt_exclude_indices", [])
        # group_map = st.session_state.get("tmt_group_map", {})
        # Get group mapping from parameters
        params = workflow_params
        group_map = {}
        for key, value in params.items():
            if key.startswith("TMT-group-") and value:
                # Extract the numeric part from keys like "TMT-group-sample1"
                match = re.search(r'sample(\d+)', key)
                if match:
                    # Subtract 1 to convert to a 0-based index (0, 1, 2...).
                    # If your samples are already 0-based, remove the -1 adjustment.
                    index = str(int(match.group(1)) - 1)
                    group_map[index] = value

        # 1. Extract keys labeled as "skip" from group_map as integer list
        exclude_indices = [
            int(k) for k, v in group_map.items() if v.lower() == "skip"
        ]

        # 2. Remove "skip" entries from group_map (keep only actual group info)
        group_map = {
            int(k): v for k, v in group_map.items() if v.lower() != "skip"
        }

        start_column_offset = 4

        # st.write("exclude_indices:", exclude_indices)
        # st.write("group_map:", group_map)

        if exclude_indices:
            # st.write("Current columns:", df.columns.tolist())
            # st.write("Number of columns:", len(df.columns))
            # st.write("Exclude indices:", exclude_indices)
            # st.write("Offset:", start_column_offset)
            cols_to_drop = [df.columns[i + start_column_offset] for i in exclude_indices]
            df_cleaned = df.drop(columns=cols_to_drop)
        else:
            df_cleaned = df.copy()

        current_cols = df_cleaned.columns.tolist()
        sample_cols = current_cols[start_column_offset:]

        # Ensure sample columns are numeric for downstream preprocessing/statistics.
        pivot_df = df_cleaned.copy()
        if sample_cols:
            pivot_df[sample_cols] = pivot_df[sample_cols].apply(pd.to_numeric, errors='coerce')

        protein_col = pivot_df.columns[0]
        expr_df = pivot_df.set_index(protein_col)[sample_cols]
        expr_df = expr_df.replace(0, np.nan)
        expr_df = np.log2(expr_df + 1)
        expr_df = expr_df.dropna()

        return pivot_df, expr_df, group_map


def get_abundance_data(workspace: Path, include_contaminants: bool = False) -> tuple | None:
    """Wrapper that handles cache key (workspace + CSV mtime + params mtime).

    Args:
        workspace: Path to the workspace directory
        include_contaminants: Keep the proteins with CON__ accessions (cRAP).
            They are left out by default, so downstream analysis never sees them.

    Returns:
        Tuple of (pivot_df, expr_df, group_map) or None if data unavailable
    """
    workflow_dir = get_workflow_dir(workspace)
    run_id = get_run_id(workspace)
    quant_dir = get_results_dir(workspace, run_id) / "quant_results"

    if not quant_dir.exists():
        return None

    csv_files = sorted(quant_dir.glob("*.csv"))
    if not csv_files:
        return None

    csv_mtime = csv_files[0].stat().st_mtime

    params_file = workflow_dir / "params.json"
    params_mtime = params_file.stat().st_mtime if params_file.exists() else 0.0

    result = load_abundance_data(str(workspace), csv_mtime, params_mtime, run_id)
    if result is None or include_contaminants:
        return result
    pivot_df, expr_df, group_map = result
    id_col = pivot_df.columns[0]
    return (
        drop_contaminants(pivot_df, id_col),
        expr_df[~expr_df.index.to_series().map(is_contaminant).to_numpy(dtype=bool)],
        group_map,
    )


def get_id_column(workspace: Path, pivot_df: pd.DataFrame) -> str:
    """Resolve the protein/row identifier column for the active analysis mode.

    LFQ reports always use "ProteinName"; TMT reports use whatever the
    report's first column is actually named (e.g. "protein").
    """
    analysis_mode = get_run_params(workspace).get("analysis-mode", "LFQ")
    return "ProteinName" if analysis_mode == "LFQ" else pivot_df.columns[0]


def get_sample_group_map(workspace: Path, pivot_df: pd.DataFrame, group_map: dict) -> dict:
    """Normalize group_map into {actual_sample_column_name: group_name}.

    LFQ group_map keys are already clean sample names (optionally with a
    ".mzML" suffix). TMT group_map keys are 0-based channel indices that must
    be matched against the report's actual "sampleN[...]" column names.
    """
    analysis_mode = get_run_params(workspace).get("analysis-mode", "LFQ")

    if analysis_mode == "LFQ":
        grouped = {
            k[:-5] if k.endswith(".mzML") else k: v
            for k, v in group_map.items()
        }
        # Every sample column is listed, those without a group as "Unassigned"
        return {
            c: grouped.get(c) or UNASSIGNED
            for c in pivot_df.columns
            if c not in ("ProteinName", "PeptideSequence", *STAT_COLUMNS)
        }

    actual_sample_names = pivot_df.columns.tolist()
    norm_map = {}
    for k, v in group_map.items():
        try:
            sample_idx = int(k) + 1
        except (TypeError, ValueError):
            continue
        target_substring = f"sample{sample_idx}["
        real_full_name = next((name for name in actual_sample_names if target_substring in name), None)
        if real_full_name:
            norm_map[real_full_name] = v if v and v.strip() else UNASSIGNED
    # Channels left without a group (not skipped) count as "Unassigned"
    for name in actual_sample_names:
        if re.search(r"sample\d+\[", str(name)) and name not in norm_map:
            norm_map[name] = UNASSIGNED
    return norm_map


def sample_columns(df: pd.DataFrame, sample_group_map: dict) -> list:
    """The sample columns of ``df``, in table order (no id or summary columns)."""
    return [c for c in df.columns if c in sample_group_map]


def group_metadata(sample_group_map: dict, columns) -> pl.DataFrame:
    """Sample-to-group table for the statistics, leaving out unassigned samples."""
    return pl.DataFrame(
        [
            {"sample_id": s, "group": sample_group_map[s]}
            for s in columns
            if sample_group_map.get(s, UNASSIGNED) != UNASSIGNED
        ],
        schema={"sample_id": pl.String, "group": pl.String},
    )


# Session-state keys of the postprocessing chain, in pipeline order.
POSTPROCESSING_STEPS = ["filtered_df", "imputed_df", "normalized_df", "statistics_df"]


def clear_downstream_steps(step: str) -> None:
    """Drop the outputs of every postprocessing step after ``step``.

    Re-applying a step changes the table the later steps were computed from,
    so their stored results no longer match and must be recomputed.
    """
    # The note of how the comparison was computed no longer applies
    st.session_state.pop("comparison_settings", None)
    for key in POSTPROCESSING_STEPS[POSTPROCESSING_STEPS.index(step) + 1:]:
        st.session_state.pop(key, None)
        if key == "normalized_df":
            st.session_state.pop("normalized_transform", None)


def postprocessing_param(params: dict, key: str, options: list | None = None):
    """Return the stored value of a postprocessing widget, seeding the default.

    Workspaces created before a key was added to ``default-parameters.json``
    lack it in their ``params.json``; fall back to the shipped default and add
    it to ``params`` so ``save_params`` persists the user's choice. When
    ``options`` is given, a stored value that is no longer a valid option
    (e.g. a two-group test after a third group was added) is replaced by the
    first option.
    """
    if key not in params:
        with open("default-parameters.json", "r", encoding="utf-8") as f:
            params[key] = json.load(f)[key]
    if options is not None:
        if params[key] not in options:
            params[key] = options[0]
        if key in st.session_state and st.session_state[key] not in options:
            del st.session_state[key]
    return params[key]


_STEP_LABELS = {
    "filtered_df": "Filtering",
    "imputed_df": "Imputation",
    "normalized_df": "Normalization",
    "statistics_df": "Statistics",
}


def get_active_table(pivot_df: pd.DataFrame) -> tuple[pd.DataFrame, str, bool]:
    """Return the most-processed protein table of the downstream analysis.

    Returns ``(table, source, is_log2)``: the output of the latest of
    Filtering, Imputation and Normalization that has been applied (falling
    back to the workflow's abundance table), the name of the step it came
    from, and whether its intensities are already log2-scaled.
    """
    for key in ["normalized_df", "imputed_df", "filtered_df"]:
        df = st.session_state.get(key)
        if df is not None:
            is_log2 = (
                key == "normalized_df"
                and st.session_state.get("normalized_transform", "log2") == "log2"
            )
            return df, _STEP_LABELS[key], is_log2
    return pivot_df, "Workflow abundance table", False


def show_pipeline_banner(uses: str = "table") -> None:
    """Tell the user which downstream-analysis output a plot page is showing.

    ``uses="table"`` reports the protein table from :func:`get_active_table`;
    ``uses="statistics"`` reports the Statistics step output.
    """
    done = [label for key, label in _STEP_LABELS.items() if st.session_state.get(key) is not None]
    steps = " → ".join(
        f"{'✅' if label in done else '⬜'} {label}" for label in _STEP_LABELS.values()
    )
    if uses == "statistics":
        source = "the **Statistics** step" if "Statistics" in done else "no statistics yet"
    else:
        applied = [label for label in done if label != "Statistics"]
        source = f"the **{applied[-1]}** step output" if applied else "the **unprocessed** workflow abundance table"
    st.info(f"Showing {source}.  \nDownstream analysis: {steps}")


def log2_matrix(df: pd.DataFrame, id_col: str, sample_cols: list, is_log2: bool) -> pd.DataFrame:
    """Protein x sample matrix on a log2 scale, zeros treated as missing."""
    mat = df.set_index(id_col)[sample_cols].apply(pd.to_numeric, errors="coerce")
    if not is_log2:
        mat = np.log2(mat.where(mat > 0))
    return mat
