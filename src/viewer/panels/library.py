"""results_library as a Viewer component."""

import streamlit as st
from src.common.results_helpers import get_run_dir
from src.viewer.scope import scoped


def render(ctx) -> None:
    run_dir = get_run_dir(ctx.workspace)
    if not run_dir.exists():
        st.warning("No workspace found. Please run a workflow first.")
        return
    library_dir = run_dir / "results" / "library"
    library_tsv = library_dir / "spectral_library.tsv"
    if library_tsv.exists():
        st.success("Spectral library available for download.")

        with open(library_tsv, "rb") as f:
            st.download_button(
                label="Download Spectral Library (TSV)",
                data=f,
                file_name="spectral_library.tsv",
                mime="text/tab-separated-values",
                use_container_width=True,
                type="primary",
                key=scoped("library_download"),
            )

        # Preview option
        with st.expander("Preview Library Contents"):
            import pandas as pd

            try:
                df = pd.read_csv(library_tsv, sep="\t", nrows=100)
                st.dataframe(df, use_container_width=True)
                st.caption("Showing first 100 rows")
            except Exception as e:
                st.error(f"Could not preview: {e}")
    else:
        st.info("No spectral library in this run. Switch on **Spectral library** on **Describe experiment** to build one.")
