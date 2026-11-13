"""Identified proteins, peptides and spectra in one component (kept for saved views)."""
import streamlit as st

from src.viewer import protein_data
from src.viewer.scope import needs_results


def render(ctx) -> None:
    data = protein_data.load(ctx)
    if data is None:
        return needs_results()
    protein_data.protein_table(data)
    if data.proteins.empty:
        return
    st.markdown("**Peptides**")
    protein_data.protein_peptides(data)
    st.markdown("**PSMs**")
    protein_data.protein_psms(data)
    st.markdown("**Spectrum**")
    protein_data.psm_spectrum(data)
