"""Single run: the MS1 peaks by retention time and m/z; a click picks the spectrum."""
from openms_insight import Heatmap

from src.viewer import run_files
from src.viewer.scope import heading, scoped
from src.viz.ms_data import raw_heatmap
from src.viz.plumbing import component, only_accepted


def render(ctx) -> None:
    run = run_files.shown(ctx)
    if run is None:
        return
    head = heading()
    if not run.has_ms1():
        # An empty plot would only clear the spectrum picked in the scan table
        run.no_ms1(head.get("title", "Peak map (MS1)"), head.get("help"))
        return
    settings = raw_heatmap.settings(title=head.get("title", "Peak map (MS1)"), spectrum="scan",
                                    zoom_identifier=f"file_peak_map_zoom_{run.stem}")
    plot = component(Heatmap, run.cache_dir, run.cache_id("ms1"),
                     lambda: {"data": run.tables()["ms1"].lazy()},
                     **only_accepted(Heatmap, **settings))
    plot(state_manager=run_files.state(ctx), height=420, key=scoped(f"file_peak_map_{run.stem}"), help=head.get("help"))
