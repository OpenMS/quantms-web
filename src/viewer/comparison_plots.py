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


# The plots, drawn the same way in the Viewer and in Edit comparison


def draw_missingness(table, samples: list[str], name: str) -> None:
    from src.common.insight_tables import show_series
    from src.common.postprocessing_plots import MISSING_COLORS, MISSING_SERIES, missingness_vs_intensity

    data, _ = missingness_vs_intensity(table, samples)
    show_series(
        data, name, x_column="x", y_column="y", series_column="series", series_names=MISSING_SERIES,
        series_colors=MISSING_COLORS, x_label="Mean log2 abundance", y_label="Density",
        y_range_mode="nonnegative",
    )


def draw_imputed(before, after, samples: list[str], name: str) -> None:
    from src.common.insight_tables import show_series
    from src.common.postprocessing_plots import IMPUTED_COLORS, IMPUTED_SERIES, imputed_value_preview

    show_series(
        imputed_value_preview(before, after, samples), name, x_column="x", y_column="y", series_column="series",
        series_names=IMPUTED_SERIES, series_colors=IMPUTED_COLORS, x_label="log2 abundance", y_label="Values",
        y_range_mode="nonnegative",
    )


def draw_distributions(before, after, samples: list[str], after_is_log: bool, name: str) -> None:
    from src.common.insight_tables import show_series
    from src.common.postprocessing_plots import STAGE_COLORS, STAGE_SERIES, sample_distributions

    from src.common.samples import short_sample_names

    data = sample_distributions(before, after, samples, after_is_log)
    data["Sample"] = data["Sample"].map(short_sample_names(samples))
    show_series(
        data, name, height=max(360, 44 * len(samples) + 140),
        x_column="Sample", y_column="Intensity", series_column="Stage", series_names=STAGE_SERIES,
        series_colors=STAGE_COLORS, mode="box", orientation="h", x_label="Sample", y_label="log2 abundance",
        y_range_mode="normal",
    )


def draw_pvalues(statistics_df, name: str) -> None:
    from src.common.insight_tables import show_series
    from src.common.postprocessing_plots import pvalue_histogram

    data, _ = pvalue_histogram(statistics_df)
    show_series(
        data, name, x_column="p-value", y_column="Proteins", mode="bars", x_label="p-value", y_label="Proteins",
        show_legend=False, y_range_mode="nonnegative",
    )
