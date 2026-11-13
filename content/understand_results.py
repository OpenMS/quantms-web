"""Understand your Results: what each view and component of Explore results shows and how to read it."""
import streamlit as st

from src.common.common import page_setup
from src.viewer.catalog import LABELS, PRESETS
from src.viewer.teaching import COMPONENTS, VIEWS, anchor

page_setup()
st.title("📖 Understand your Results")
st.markdown(
    "What **Explore results** shows, view by view. Each component's title there links to its section here."
)

shown = set()
for preset in PRESETS:
    st.header(preset["label"], anchor=anchor(preset["name"]))
    st.markdown(VIEWS.get(preset["name"], preset.get("help", "")))
    for row in preset["layout"]:
        for name in row:
            if name in shown or name not in COMPONENTS:
                continue
            shown.add(name)
            st.subheader(LABELS.get(name, name), anchor=anchor(name))
            st.markdown(COMPONENTS[name])

others = [name for name in LABELS if name not in shown and name in COMPONENTS]
if others:
    st.header("More components", anchor="more-components")
    st.markdown("Components no ready-made view shows; add them to a view of your own with ✏️ or **+**.")
    for name in others:
        st.subheader(LABELS[name], anchor=anchor(name))
        st.markdown(COMPONENTS[name])
