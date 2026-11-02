"""SequenceView renders with the installed openms-insight release."""
from pathlib import Path

import polars as pl
from streamlit.testing.v1 import AppTest

PAGE = '''
import sys
sys.path.insert(0, {root!r})
import src.common.insight_compat  # noqa: F401
import polars as pl
import streamlit as st
from openms_insight import SequenceView, StateManager

sv = SequenceView(
    cache_id="compat_sv",
    sequence_data=pl.LazyFrame({{"scan_id": [1], "sequence": ["PEPTIDER"], "precursor_charge": [2]}}),
    peaks_data=pl.LazyFrame({{"scan_id": [1], "peak_id": [1], "mass": [200.5], "intensity": [1.0]}}),
    cache_path={cache!r},
    filters={{"spectrum": "scan_id"}},
)
state = StateManager("compat")
state.set_selection("spectrum", 1)
sv(key="compat_sv", state_manager=state)
st.write("rendered")
'''


def test_sequenceview_renders(tmp_path):
    page = tmp_path / "page.py"
    page.write_text(PAGE.format(root=str(Path.cwd()), cache=str(tmp_path / "cache")))
    at = AppTest.from_file(str(page), default_timeout=60).run()
    assert not at.exception, at.exception[0].value
    assert [m.value for m in at.markdown] == ["rendered"]
