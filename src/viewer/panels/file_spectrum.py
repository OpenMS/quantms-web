"""Single run: the peaks of the picked spectrum."""
import streamlit as st
from openms_insight import DeconvolvedSpectrum

from src.common.plot_style import spectrum_styling
from src.viewer import run_files
from src.viewer.scope import heading, scoped, zoom
from src.viz.ms_data import spectrum
from src.viz.plumbing import component, only_accepted


def render(ctx) -> None:
    run = run_files.shown(ctx)
    if run is None:
        return
    state = run_files.state(ctx)
    head = heading()
    if state.get_selection("scan") is None:
        st.markdown(f"**{head.get('title', 'Spectrum')}**", help=head.get("help"))
        st.info("Pick a spectrum in the scan table, the chromatogram or the peak map.")
        return
    # The library's raw spectrum, without masses: the peaks of the spectrum only
    settings = spectrum.settings("scan", styling=spectrum_styling(), signal_column=None,
                                 title=head.get("title", "Spectrum"))
    plot = component(DeconvolvedSpectrum, run.cache_dir, run.cache_id("spectra"),
                     lambda: {"data": run_files.spectra(run.tables()["peaks"]).lazy()},
                     **only_accepted(DeconvolvedSpectrum, **settings), **zoom(DeconvolvedSpectrum, "file_spectrum"))
    plot(state_manager=state, height=360, key=scoped(f"file_spectrum_{run.stem}"), help=head.get("help"))
