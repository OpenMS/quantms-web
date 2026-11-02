"""Stepping through a PSM table, and where a PSM was reached from.

Previous/Next follow the table as the viewer sees it: openms-insight's
``Table.export_data`` applies the linked filters, range and toggle filters,
the column filters and the sort the viewer chose, so the neighbours are the
rows above and below in that order. A step sets the table's selections as a
click on that row would, and the table pages to it.
"""

from __future__ import annotations

import html

import streamlit as st


def _same(a, b) -> bool:
    """Selection values compare loosely: the browser sends 3.0 for 3."""
    if a is None or b is None:
        return False
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return str(a) == str(b)


def row_position(table, state: dict, column: str, current) -> tuple[int | None, list]:
    """(index of ``current`` in ``column``, all values of ``column``), in the table's
    current filter and sort order. The index is None when the row is not shown."""
    values = table.export_data(state)[column].to_list()
    for i, value in enumerate(values):
        if _same(value, current):
            return i, values
    return None, values


def neighbour(table, state: dict, interactivity: dict, primary: str, step: int) -> dict | None:
    """Selections for the row ``step`` rows away from the selected one, or None at an end."""
    column = interactivity[primary]
    rows = table.export_data(state)
    values = rows[column].to_list()
    current = state.get(primary)
    index = next((i for i, v in enumerate(values) if _same(v, current)), None)
    if index is None or not 0 <= index + step < len(values):
        return None
    row = rows.row(index + step, named=True)
    return {identifier: row[col] for identifier, col in interactivity.items() if col in row}


def psm_stepper(table, state_manager, primary: str, key: str, interactivity: dict | None = None,
                noun: str = "PSM") -> None:
    """◀ Previous · "PSM i of n" · Next ▶ for the table's selected row.

    ``primary`` is the selection identifier that names the row (e.g.
    "identification"); ``interactivity`` maps identifiers to columns and
    defaults to the table's own.
    """
    interactivity = dict(interactivity or table._interactivity)
    state = state_manager.get_all_selections()
    index, values = row_position(table, state, interactivity[primary], state.get(primary))

    def go(step: int) -> None:
        target = neighbour(table, state_manager.get_all_selections(), interactivity, primary, step)
        for identifier, value in (target or {}).items():
            state_manager.set_selection(identifier, value)

    c1, c2, c3 = st.columns([1, 2, 1], vertical_alignment="center")
    c1.button(
        "◀ Previous", key=f"{key}-prev", on_click=go, args=(-1,),
        disabled=index is None or index == 0, use_container_width=True,
        help=f"The {noun} one row up in the table, in its current sort and filters",
    )
    if index is None:
        c2.markdown(
            f"<div style='text-align:center'>Selected {noun} not in the table's current view "
            f"({len(values):,} rows)</div>", unsafe_allow_html=True,
        )
    else:
        c2.markdown(
            f"<div style='text-align:center'>{noun} <b>{index + 1:,}</b> of {len(values):,} "
            "in the table's order</div>", unsafe_allow_html=True,
        )
    c3.button(
        "Next ▶", key=f"{key}-next", on_click=go, args=(1,),
        disabled=index is None or index + 1 >= len(values), use_container_width=True,
        help=f"The {noun} one row down in the table, in its current sort and filters",
    )


def breadcrumb(parts: list[tuple[str, str | None]]) -> None:
    """``Proteins › P00918 › PEPTIDEK › scan 1234``: each part (label, anchor) links to
    its section on the page; the last one, the current item, is plain text."""
    crumbs = []
    for i, (label, anchor) in enumerate(parts):
        text = html.escape(label)
        if anchor and i < len(parts) - 1:
            crumbs.append(f"<a href='#{anchor}' target='_self'>{text}</a>")
        else:
            crumbs.append(f"<b>{text}</b>")
    st.markdown(
        "<nav aria-label='breadcrumb' style='font-size:0.95rem;margin:0.25rem 0 0.5rem'>"
        + " <span style='opacity:0.6'>›</span> ".join(crumbs)
        + "</nav>",
        unsafe_allow_html=True,
    )
