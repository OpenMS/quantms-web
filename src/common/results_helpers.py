"""Helper functions for results pages."""
import re
import pandas as pd
import polars as pl
from pathlib import Path
from pyopenms import IdXMLFile, PeptideIdentificationList, MSExperiment, MzMLFile


def get_workflow_dir(workspace):
    """Get the workflow directory path."""
    return Path(workspace, "topp-workflow")


def idxml_to_df(idxml_file):
    """Parse idXML file and return DataFrame with peptide hits."""
    proteins = []
    peptides = PeptideIdentificationList()
    IdXMLFile().load(str(idxml_file), proteins, peptides)
    peptides = [peptides.at(i) for i in range(peptides.size())]

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


def extract_filename_from_idxml(idxml_path: Path) -> str:
    """Derive mzML filename from idXML filename."""
    stem = idxml_path.stem
    for suffix in ['_comet', '_per', '_filter']:
        stem = stem.replace(suffix, '')
    return f"{stem}.mzML"


def parse_idxml(idxml_path: Path) -> tuple[pl.DataFrame, list[str]]:
    """Parse idXML and return DataFrame for openms_insight.

    Returns:
        Tuple of (id_df, spectra_data list of source filenames)
    """
    proteins = []
    peptides = PeptideIdentificationList()
    IdXMLFile().load(str(idxml_path), proteins, peptides)
    peptides = [peptides.at(i) for i in range(peptides.size())]

    # Derive mzML filename from idXML filename (e.g., 02COVID_filter.idXML -> 02COVID.mzML)
    spectra_data = [extract_filename_from_idxml(idxml_path)]

    # Build filename to index mapping
    filename_to_index = {Path(f).name: i for i, f in enumerate(spectra_data)}

    records = []
    for pep in peptides:
        # Get spectrum reference from meta value (key may be bytes or string)
        spec_ref = ""
        if pep.metaValueExists("spectrum_reference"):
            spec_ref = pep.getMetaValue("spectrum_reference")
            if isinstance(spec_ref, bytes):
                spec_ref = spec_ref.decode()
        scan_id = extract_scan_from_ref(spec_ref)

        # Get file index from id_merge_index or derive from filename
        file_index = pep.getMetaValue("id_merge_index") if pep.metaValueExists("id_merge_index") else 0
        filename = spectra_data[file_index] if file_index < len(spectra_data) else ""

        for h in pep.getHits():
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
                "protein_accession": ";".join([ev.getProteinAccession() for ev in h.getPeptideEvidences()]),
            })

    return pl.DataFrame(records), spectra_data


def build_spectra_cache(mzml_dir: Path, filename_to_index: dict) -> tuple[pl.DataFrame, dict]:
    """Extract MS2 spectra from mzML files and return DataFrame.

    Args:
        mzml_dir: Directory containing mzML files
        filename_to_index: Dict mapping filename to file_index

    Returns:
        Tuple of (spectra_df, updated filename_to_index)
    """
    records = []
    peak_id = 0

    for mzml_path in sorted(mzml_dir.glob("*.mzML")):
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
