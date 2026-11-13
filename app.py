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
    # One section follows the path a user takes: describe the experiment, run
    # it, then read every result on Explore results.
    online = st.session_state.settings["online_deployment"]
    export, export_icon = ("Download results", "⬇️") if online else ("Export results", "📤")
    pages = {
        "OpenDDA": [
            st.Page(Path("content", "quickstart.py"), title="Quickstart", icon="👋", default=True),
        ],
        "🧬 DDA analysis": [
            st.Page(Path("content", "workflow_configure.py"), title="Describe experiment", icon="⚙️"),
            st.Page(Path("content", "workflow_run.py"), title="Run analysis", icon="🚀"),
            st.Page(Path("content", "viewer.py"), title="Explore results", icon="👀"),
            # A local install exports results to a folder instead of downloading them
            st.Page(Path("content", "workflow_export.py"), title=export, icon=export_icon),
            st.Page(Path("content", "understand_results.py"), title="Understand your Results", icon="📖"),
        ],
    }

    pg = st.navigation(pages)
    pg.run()

