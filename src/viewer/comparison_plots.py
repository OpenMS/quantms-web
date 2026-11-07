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
