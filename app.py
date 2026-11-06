import os
# Polars' default (CPU-core-count-sized) native thread pool access-violation
# crashes the whole process on some high-core-count Windows machines when
# invoked from Streamlit's script-runner thread. Must be set before polars
# is imported anywhere (including transitively via openms_insight).
os.environ.setdefault("POLARS_MAX_THREADS", "1")

import streamlit as st
from pathlib import Path
import json
# For some reason the windows version only works if this is imported here
import pyopenms
# Patches openms-insight classes before any page renders them
import src.common.insight_compat  # noqa: F401

if "settings" not in st.session_state:
        with open("settings.json", "r") as f:
            st.session_state.settings = json.load(f)

if __name__ == '__main__':
    # Sections follow the path a user takes: set up and run, then read every
    # result in the Viewer.
    pages = {
        "Start here": [
            st.Page(Path("content", "quickstart.py"), title="Home", icon="🏠", default=True),
        ],
        "Step 1 · Set up and run": [
            st.Page(Path("content", "workflow_configure.py"), title="Describe experiment", icon="⚙️"),
            st.Page(Path("content", "workflow_run.py"), title="Run analysis", icon="🚀"),
        ],
        "Step 2 · Results": [
            st.Page(Path("content", "viewer.py"), title="Viewer", icon="👁️"),
        ],
    }

    pg = st.navigation(pages)
    pg.run()

