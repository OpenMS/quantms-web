"""Single run: the scan table of the run; a click picks the spectrum."""
from openms_insight import Table

from src.viewer import run_files
from src.viewer.scope import heading, scoped
from src.viz.ms_data import scan_table
from src.viz.plumbing import component


def render(ctx) -> None:
    run = run_files.shown(ctx)
    if run is None:
        return
    head = heading()
    table = component(Table, run.cache_dir, run.cache_id("scans"),
                      lambda: {"data": run.tables()["scans"].lazy()},
                      **scan_table.settings("scan", clears=(), peaks=True, title=head.get("title", "Scans")))
    table(state_manager=run_files.state(ctx), height=376, key=scoped(f"file_scans_{run.stem}"), help=head.get("help"))
    run.note()
