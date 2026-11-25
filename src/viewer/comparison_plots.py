"""What the plots behind the steps of the group comparison draw from.

Each plot is a Viewer component of its own (src/viewer/panels/comparison_*.py).

They show the comparison as computed (recommended settings, or what was applied
in Edit comparison), so they are read alongside the volcano plot and the table.
Edit comparison shows the same plots as previews of a setting before it is
applied.
"""

import streamlit as st

from src.common.guide import ensure_group_comparison
from src.common.results_helpers import (
    get_abundance_data, get_sample_group_map, group_metadata, sample_columns,
)
from src.viewer.scope import needs_results


def comparison(ctx):
    """(step tables, compared samples, sample -> group), or None after a note."""
    result = get_abundance_data(ctx.workspace)
    if result is None:
        needs_results()
        return None
    if ensure_group_comparison(ctx.params) is None:
        return None
    pivot_df, _, group_map = result
    sample_group_map = get_sample_group_map(ctx.workspace, pivot_df, group_map)
    filtered = st.session_state.get("filtered_df")
    columns = sample_columns(filtered if filtered is not None else pivot_df, sample_group_map)
    samples = group_metadata(sample_group_map, columns)["sample_id"].to_list()
    return st.session_state, samples, sample_group_map


# The plots, drawn the same way in the Viewer and in Edit comparison: the
# library's panels (src/viz/comparison_steps.py), in this app's cache and state


def plot_height(samples: list[str]) -> int:
    from src.viz.comparison_steps import plot_height as height

    return height(samples)


def _where(name: str) -> dict:
    """Cache, widget key, and the title and help the Viewer gives the component drawn."""
    from src.common.insight_tables import _default_cache_dir
    from src.viewer.scope import heading, scoped

    head = heading()
    return dict(cache_path=_default_cache_dir(), key=scoped(f"series_{name}"), title=head.get("title") or "",
                help=head.get("help"))


def draw_missingness(table, samples: list[str], name: str) -> None:
    from src.viz.comparison_steps import show_missing_values

    show_missing_values(table, samples, name, **_where(name))


def draw_imputed(before, after, samples: list[str], name: str) -> None:
    from src.viz.comparison_steps import show_imputed_values

    show_imputed_values(before, after, samples, name, **_where(name))


def draw_distributions(before, after, samples: list[str], after_is_log: bool, name: str) -> None:
    from src.common.samples import short_sample_names
    from src.viz.comparison_steps import show_normalization

    show_normalization(before, after, samples, after_is_log, name, labels=short_sample_names(samples),
                       **_where(name))


def draw_pvalues(statistics_df, name: str, height: int = 360) -> None:
    from src.viz.comparison_steps import show_pvalues

    show_pvalues(statistics_df, name, noun="Proteins", height=height, **_where(name))
