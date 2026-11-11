"""Accept or reject the PSM shown, with a note, and the list of verdicts of a run."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import streamlit as st

from src.common.psm_curation import load_verdicts, set_verdict, verdict_of, verdicts_tsv
from src.common.spectrum_view import psm_lookup

_LABELS = {"accepted": "Accept", "rejected": "Reject"}
_VERDICT = {label: verdict for verdict, label in _LABELS.items()}


def _hits(idxml: Path) -> pd.DataFrame:
    return psm_lookup(str(idxml), idxml.stat().st_mtime)


def _psm(hit) -> dict:
    return {"file": hit["filename"], "scan": int(hit["scan_id"]), "sequence": hit["sequence"],
            "charge": int(hit["charge"]), "mz": float(hit["mz"]), "rt": float(hit["rt"])}


def _widget_key(key: str, psm: dict) -> str:
    name = f"{psm['file']}|{psm['scan']}|{psm['sequence']}|{psm['charge']}"
    return f"{key}_{hashlib.sha1(name.encode()).hexdigest()[:12]}"


def _save(results_dir: Path, psm: dict, widget: str) -> None:
    verdict = _VERDICT.get(st.session_state.get(f"{widget}_verdict"))
    note = st.session_state.get(f"{widget}_note", "")
    set_verdict(results_dir, psm, verdict, note)


def verdict_controls(results_dir: Path, idxml: Path, state_manager, key: str) -> None:
    """Accept / Reject and a note for the PSM selected in the run's PSM table."""
    psm_id = state_manager.get_selection("identification")
    hits = _hits(idxml)
    if psm_id is None or int(psm_id) not in hits.index:
        return
    psm = _psm(hits.loc[int(psm_id)])
    verdict, note = verdict_of(load_verdicts(results_dir), psm)
    widget = _widget_key(key, psm)
    left, right = st.columns([2, 5], vertical_alignment="bottom")
    with left:
        st.segmented_control(
            "Your verdict on this PSM", list(_VERDICT), default=_LABELS.get(verdict), key=f"{widget}_verdict",
            on_change=_save, args=(results_dir, psm, widget),
            help="Kept with the run (results/curation.tsv) and in its ZIP. Click the chosen verdict again to clear it.",
        )
    with right:
        st.text_input(
            "Note", value=note, key=f"{widget}_note", on_change=_save, args=(results_dir, psm, widget),
            placeholder="Why, for example: few fragments matched, an unexplained intense peak",
        )


def verdict_list(results_dir: Path, idxml: Path, state_manager, key: str, run_name: str = "run") -> None:
    """The run's verdicts, filtered by verdict, with a TSV download; a row picks its PSM."""
    frame = load_verdicts(results_dir)
    accepted = int((frame["verdict"] == "accepted").sum())
    rejected = int((frame["verdict"] == "rejected").sum())
    label = f"Your verdicts ({accepted} accepted, {rejected} rejected)" if len(frame) else "Your verdicts"
    with st.expander(label, expanded=False):
        if frame.empty:
            st.caption(
                "Pick a PSM in the table and accept or reject it under its spectrum. "
                "Your verdicts and notes are listed here and can be downloaded as a table."
            )
            return
        show = st.segmented_control("Show", ["All", "Accepted", "Rejected"], default="All", key=f"{key}_show")
        if show in ("Accepted", "Rejected"):
            frame = frame[frame["verdict"] == show.lower()]
        shown = frame.assign(Verdict=frame["verdict"].str.capitalize())
        event = st.dataframe(
            shown[["Verdict", "sequence", "charge", "note", "scan", "file", "mz", "rt"]].rename(columns={
                "sequence": "Peptide", "charge": "Charge", "file": "File", "scan": "Scan", "mz": "m/z",
                "rt": "RT (s)", "note": "Note"}),
            hide_index=True, width="stretch", on_select="rerun", selection_mode="single-row",
            key=f"{key}_list_{st.session_state.get(f'{key}_picks', 0)}", height=min(38 + 35 * len(shown), 400),
            column_config={"Peptide": st.column_config.TextColumn(width="medium"),
                           "Note": st.column_config.TextColumn(width="large"),
                           "File": st.column_config.TextColumn(width="small"),
                           "m/z": st.column_config.NumberColumn(format="%.4f"),
                           "RT (s)": st.column_config.NumberColumn(format="%.1f")},
        )
        rows = event.selection.rows if event is not None else []
        if rows:
            _pick(shown.iloc[rows[0]], idxml, state_manager, key)
        st.download_button(
            "Download verdicts (TSV)", verdicts_tsv(load_verdicts(results_dir)),
            file_name=f"{run_name}-psm-verdicts.tsv", mime="text/tab-separated-values", key=f"{key}_tsv",
            icon=":material/download:",
        )


def _pick(row, idxml: Path, state_manager, key: str) -> None:
    """Show the PSM of a list row in the table, the spectrum and the views linked to them."""
    hits = _hits(idxml)
    match = hits[(hits["filename"] == row["file"]) & (hits["scan_id"] == int(row["scan"]))
                 & (hits["sequence"] == row["sequence"]) & (hits["charge"] == int(row["charge"]))]
    if match.empty:
        st.caption(f"{row['sequence']} (scan {int(row['scan'])}) is not among this file's PSMs. "
                   "Pick its file above to show it.")
        return
    psm_id = int(match.index[0])
    state_manager.set_selection("identification", psm_id)
    state_manager.set_selection("file", int(match.iloc[0]["file_index"]))
    state_manager.set_selection("spectrum", int(row["scan"]))
    # A new list, its row unticked, so that filtering it later never picks another PSM
    st.session_state[f"{key}_picks"] = st.session_state.get(f"{key}_picks", 0) + 1
    st.rerun()
