"""Single run: total ion current and base peak of the MS1 spectra; a click picks the spectrum."""
from openms_insight import SeriesPlot

from src.viewer import run_files
from src.viewer.scope import heading, scoped, zoom
from src.viz.ms_data import chromatogram
from src.viz.plumbing import component, only_accepted

SERIES = ["TIC", "Base peak"]


def render(ctx) -> None:
    run = run_files.shown(ctx)
    if run is None:
        return
    head = heading()
    settings = chromatogram.settings(SERIES, level=1, spectrum="scan")
    if head.get("title"):
        settings["title"] = head["title"]
    plot = component(SeriesPlot, run.cache_dir, run.cache_id("chromatogram"),
                     lambda: {"data": run.tables()["chromatogram"].lazy()},
                     **only_accepted(SeriesPlot, **settings), **zoom(SeriesPlot, "file_chromatogram"))
    plot(state_manager=run_files.state(ctx), height=320, key=scoped(f"file_chromatogram_{run.stem}"),
         help=head.get("help"))
