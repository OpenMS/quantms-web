"""Identification diagnostics as a Viewer component: every file's ID rate, the files
that fall short first, with their likely causes."""

from src.common.id_diagnosis import diagnose_identifications
from src.common.insight_tables import show_insight_table
from src.viewer.panels.overview import load_overview
from src.viewer.scope import needs_results


def render(ctx) -> None:
    loaded = load_overview(ctx)
    if loaded is None:
        return needs_results()
    per_file, _, workflow_params, cache_dir, sources = loaded
    table = diagnose_identifications(per_file, workflow_params, all_files=True)
    for col in ("ID rate (%)", "Expected (%)"):
        table[col] = table[col].map(lambda v: "" if v is None or v != v else f"{v:g}")
    show_insight_table(
        table,
        name="id_diagnostics",
        cache_dir=cache_dir,
        source_files=sources,
        # The expected rate follows the search settings, which the sources do not show
        extra_key=str(table["Expected (%)"].iloc[0]) if len(table) else "",
        title="Identification diagnostics",
        # A row with causes wraps to about three lines
        height=min(640, 120 + sum(110 if isinstance(c, str) and c.strip() else 40 for c in table["Likely causes"])),
        layout="fitColumns",
        # The causes are sentences: wrap them in a wide column rather than cut them off
        column_overrides={
            "File": {"widthGrow": 2},
            "ID rate (%)": {"widthGrow": 1},
            "Expected (%)": {"widthGrow": 1},
            "Likely causes": {"formatter": "textarea", "widthGrow": 7},
        },
    )
