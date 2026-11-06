import streamlit as st
from pathlib import Path
import re
import json
import shutil
import sys
import pandas as pd
import plotly.express as px
from streamlit_plotly_events import plotly_events
from scipy.stats import ttest_ind
import numpy as np
import mygene

from collections import defaultdict        
from scipy.stats import fisher_exact         
from src.workflow.WorkflowManager import WorkflowManager
from src.common.common import page_setup
from src.common.results_helpers import get_abundance_data
from src.common.results_helpers import load_idxml, search_results_dir, write_spectra_cache
from src.common.identified_spectra import SPECTRA_DIR, write_identified_spectra
from src.common.psm_views import build_psm_views
from src.common.contaminants import CONTAMINANT_PREFIX, append_fasta, write_contaminants
from src.common.fdr_display import FDR_QUICK_CHOICES, fdr_applied, fdr_label, matching_quick_choice, percent_label, q_label
from src.common.isobaric import ISOBARIC_CHANNELS
from src.common.ms2_prediction import peptdeep_available
from src.common.match_between_runs import MATCH_BETWEEN_RUNS_KEY, lfq_transfer_params, match_between_runs
from src.common.lfq_design import groups_from_design, initial_groups, is_fractionated, seed_design
from src.common.qpx_export import app_version, export_run
from src.workflow.ExperimentalDesign import design_path, load_design, save_design, to_openms_design
from src.common import instrument
from src.pride_ui import pride_import_section

# Accession prefixes that mark decoy entries in a FASTA database.
DECOY_PREFIXES = ("DECOY_", "decoy_", "rev_", "REV_", "XXX_", "reversed_")


@st.cache_data(show_spinner=False)
def detect_decoy_prefix(fasta_path: str, mtime: float) -> str | None:
    """Return the decoy prefix used in a FASTA database, or None if it has no decoys.

    `mtime` is only part of the cache key, so a replaced file is scanned again.
    """
    try:
        with open(fasta_path, "r", errors="ignore") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:]
                    for prefix in DECOY_PREFIXES:
                        if header.startswith(prefix):
                            return prefix
    except OSError:
        return None
    return None


# Plain-language help shown above the search settings (LFQ and TMT).
DECOY_INFO = """
**Why decoys?** To estimate how many matches are wrong, the search also looks for
peptides from fake ("decoy") proteins: your proteins with their sequences reversed or
shuffled. Every decoy hit is a known false match, so the number of decoy hits tells how
many false matches hide among the real ones. *Reverse* is the standard choice.
"""

SEARCH_INFO = """
**How peptides are identified.** ProSE, the OpenMS search engine, compares every
fragmentation (MS2) spectrum with the fragment ions of the peptides of your FASTA database.
The defaults fit a tryptic digest with carbamidomethylated cysteines. Pick your **Instrument**
in *Your analysis* above to set the tolerances for your mass spectrometer; change the
rest only if your sample preparation differs:

* **enzyme**: the protease you digested with (usually Trypsin/P).
* **missed_cleavages**: how many cut sites the enzyme may skip within one peptide (2 is typical).
* **fixed**: modifications present on every copy of an amino acid, e.g.
  *Carbamidomethyl (C)* after reduction and alkylation with iodoacetamide.
* **variable**: modifications that may or may not be present, e.g. *Oxidation (M)*
  or *Phospho (STY)*. Each one added makes the search slower and less sensitive, so add only
  those you expect.
* **precursor / fragment mass tolerance**: how far a measured mass may lie from the expected
  one. Set by the Instrument. A precursor tolerance above 1 Da searches openly for
  unexpected mass shifts (*Search: Open* in *Your analysis*).
"""

SEARCH_DECOY_INFO = """* **Decoys**: your FASTA file already has decoy proteins (accessions starting with
  *{prefix}*), so they are searched as they are.
"""

# ProSE settings shown by default; the rest is in the tool's "All N parameters" dialog.
SEARCH_INCLUDE = [
    "Search:enzyme", "Search:peptide:missed_cleavages",
    "Search:modifications:fixed", "Search:modifications:variable",
    "Search:precursor:mass_tolerance_lower", "Search:precursor:mass_tolerance_upper",
    "Search:precursor:mass_tolerance_unit", "Search:precursor:isotopes",
    "Search:fragment:mass_tolerance", "Search:fragment:mass_tolerance_unit",
]
# Decoys come from the database step, rescoring and FDR from the steps after the search.
SEARCH_EXCLUDE = [
    "Search:decoys", "Search:decoy_prefix", "Search:FDR:PSM", "Search:FDR:protein",
    "percolator_executable",
]
# Orbitrap HCD (the "Orbitrap HCD" preset), the most common DDA setup.
SEARCH_DEFAULTS = {
    "Search:enzyme": "Trypsin/P",
    "Search:decoys": "auto",
    "Search:peptide:missed_cleavages": 2,
    "Search:peptide:min_size": 6,
    "Search:peptide:max_size": 40,
    "Search:modifications:fixed": "Carbamidomethyl (C)",
    "Search:modifications:variable": "Oxidation (M)\nAcetyl (Protein N-term)",
    "Search:modifications:variable_max_per_peptide": 3,
    "Search:precursor:mass_tolerance_lower": 20.0,
    "Search:precursor:mass_tolerance_upper": 20.0,
    "Search:precursor:mass_tolerance_unit": "ppm",
    "Search:precursor:isotopes": "0\n1",
    "Search:precursor:min_charge": 2,
    "Search:precursor:max_charge": 4,
    "Search:fragment:mass_tolerance": 20.0,
    "Search:fragment:mass_tolerance_unit": "ppm",
    "Search:report:top_hits": 1,
    "Search:FDR:PSM": 0.0,
    "Search:FDR:protein": 0.0,
}
SEARCH_DEFAULTS_TMT = {
    **SEARCH_DEFAULTS,
    "Search:precursor:mass_tolerance_lower": 4.5,
    "Search:precursor:mass_tolerance_upper": 4.5,
    "Search:modifications:variable": "Oxidation (M)\nAcetyl (Protein N-term)\nTMT6plex (K)\nTMT6plex (N-term)",
}


PERCOLATOR_INFO = """
**Telling right from wrong matches.** Percolator learns, from many features of each
match (ProSE's score, mass errors, matched ions, peptide length, ...), what distinguishes real matches from the
decoy ones, and re-scores every match. This typically identifies 20 to 50 % more peptides
at the same error rate. The defaults suit every experiment.
"""

FDR_INFO = """
**How many wrong matches you accept.** The false discovery rate (FDR) is the share of
reported matches expected to be wrong. **1 %** (q ≤ 0.01) is the community standard and what journals
expect; raise it only for exploration. 100 % switches the filter off.
"""


def demo_group(filename: str) -> str:
    """Group of a demo file (UPS2 spiked in at two ratios), so the demo runs as shipped."""
    if not Path("example-data", "mzML", filename).is_file():
        return ""
    for ratio in ("ratio01", "ratio10"):
        if ratio in filename:
            return f"UPS2_{ratio}"
    return ""


# Inputs the runs can come as: vendor files are converted to mzML when a run starts.
MS_DATA_TYPES = ["mzML", "raw", "d"]
DATABASE_TYPES = ("fasta", "fa")
# The demo runs, each in the folder of its group, and the demo database.
EXAMPLE_MZML = {
    str(f): demo_group(f.name) for f in sorted(Path("example-data", "mzML").glob("*.mzML"))
}
EXAMPLE_FASTA = [str(f) for f in sorted(Path("example-data", "db").glob("*.fasta"))]

MODE_LABELS = {"LFQ": "Label-free (LFQ)", "TMT": "Isobaric labels (TMT, iTRAQ)"}
MODE_HELP = {
    "LFQ": "every sample was measured in its own run (or its fractions).",
    "TMT": "samples were labelled with isobaric tags and pooled into one run (or its fractions).",
}
INSTRUMENT_LABELS = {
    instrument.ORBITRAP: "Orbitrap (HCD)",
    instrument.TIMSTOF: "timsTOF",
    instrument.ION_TRAP: "Ion trap (CID)",
}
OPEN_SEARCH, CLOSED_SEARCH = "Open Search", "Closed Search"
SEARCH_LABELS = {CLOSED_SEARCH: "Standard", OPEN_SEARCH: "Open (±500 Da)"}
SEARCH_HELP = (
    "**Standard**: the precursor mass has to match a peptide of the database with the "
    "modifications you chose, within the instrument's tolerance. **Open**: any precursor "
    "within ±500 Da matches, so peptides with modifications nobody searched for are found "
    "too; see *Mass shifts* for what they are."
)


def _no_preset_buttons(*args, **kwargs):
    """The instrument and search presets are chosen in Your analysis, not by buttons.

    A module-level function: the workflow object is pickled for its run.
    """
    return None


@st.cache_data(show_spinner="Looking at the selected runs...", max_entries=64)
def _instrument_recommendation(paths: tuple, converted_dir: str, mtimes: tuple):
    """Cached per file set; `mtimes` makes a replaced file count again."""
    return instrument.recommend(list(paths), converted_dir)


# params = page_setup()
# Threads one ProSE process uses well: on a machine with many cores, several
# processes at 16-32 threads, each searching a share of the files, finish
# sooner than one process with all of them (ProSE's documentation).
SEARCH_THREADS_PER_PROCESS = 16


def search_batches(n_files: int, max_threads: int) -> list[list[int]]:
    """Split the files of a search into the ProSE processes that run side by side.

    One process per SEARCH_THREADS_PER_PROCESS threads of the budget (at
    least one, at most one per file); the files are dealt out in order, as
    evenly as possible, so each process gets a share of the threads.
    """
    if n_files <= 0:
        return []
    processes = max(1, min(n_files, int(max_threads) // SEARCH_THREADS_PER_PROCESS))
    size, extra = divmod(n_files, processes)
    batches, start = [], 0
    for p in range(processes):
        end = start + size + (1 if p < extra else 0)
        batches.append(list(range(start, end)))
        start = end
    return batches


class WorkflowTest(WorkflowManager):
    # Every run is kept as its own result set (runs/<id>/); the results pages
    # show the run picked under "Showing run".
    multi_run = True

    def __init__(self) -> None:
        super().__init__("TOPP Workflow", st.session_state["workspace"])
        self.ui.preset_buttons = _no_preset_buttons

    def upload(self) -> None:
        """The upload widgets (the Browse dialogs of Your analysis draw the same ones)."""
        t = st.tabs(["MS data", "Protein database (FASTA)"])
        with t[0]:
            self.ui.upload_widget(
                key="mzML-files",
                name="MS data",
                file_types=MS_DATA_TYPES,
                fallback=EXAMPLE_MZML,
                suggest_folders=True,
            )
        with t[1]:
            self.ui.upload_widget(
                key="fasta-file",
                name="Protein database",
                file_types=DATABASE_TYPES,
                fallback=EXAMPLE_FASTA,
            )

    def current_params(self) -> dict:
        """The parameters with the values the widgets hold in this rerun
        (self.params is stale within a fragment rerun)."""
        prefix = self.parameter_manager.param_prefix
        current = dict(self.params)
        current.update({k[len(prefix):]: v for k, v in st.session_state.items() if k.startswith(prefix)})
        return current

    def mode(self) -> str:
        """LFQ or TMT, as chosen in Your analysis."""
        mode = self.current_params().get("analysis-mode", "LFQ")
        return mode if mode in MODE_LABELS else "LFQ"

    def _choice(self, key, options, labels, default, label, help=None):
        """A segmented choice stored as the workflow parameter `key`."""
        widget = f"{self.parameter_manager.param_prefix}{key}"
        if st.session_state.get(widget) not in options:
            stored = self.params.get(key, default)
            st.session_state[widget] = stored if stored in options else default
        st.segmented_control(label, options, format_func=lambda o: labels[o], key=widget, help=help,
                             on_change=self._keep_choice, args=(widget, default))
        self.parameter_manager.save_parameters()
        return st.session_state[widget]

    @staticmethod
    def _keep_choice(widget, default):
        # A click on the chosen option clears a segmented control: keep a choice
        if st.session_state.get(widget) is None:
            st.session_state[widget] = default

    # --- Instrument and search: presets of presets.json, read back from the settings ---

    @staticmethod
    def search_instance(mode: str) -> str:
        return "ProSE-TMT" if mode == "TMT" else "ProSE"

    def _search_values(self, mode: str) -> dict:
        """The search settings in effect: defaults, the tool's ini and the stored values."""
        instance = self.search_instance(mode)
        values = dict(SEARCH_DEFAULTS_TMT if mode == "TMT" else SEARCH_DEFAULTS)
        values.update(self.parameter_manager.get_topp_parameters("ProSE", tool_instance_name=instance))
        # Fresh from params.json: a preset applied in this rerun is not in self.params
        values.update(self.parameter_manager.get_parameters_from_json().get(instance, {}))
        return values

    @staticmethod
    def _is_open(values: dict) -> bool:
        return (values.get("Search:precursor:mass_tolerance_unit") == "Da"
                and float(values.get("Search:precursor:mass_tolerance_upper", 0)) >= 100)

    def current_instrument(self, mode: str) -> str | None:
        """The instrument preset whose tolerances are set, or None for tolerances of your own."""
        values = self._search_values(mode)
        is_open = self._is_open(values)
        same = self.parameter_manager._same_value
        presets = self.parameter_manager.load_presets()
        for name in instrument.PRESETS:
            preset = presets.get(name, {}).get(self.search_instance(mode), {})
            if not preset:
                continue
            if all(
                same(values.get(k), v) for k, v in preset.items()
                if not (is_open and k.startswith("Search:precursor:"))
            ):
                return name
        return None

    def current_search(self, mode: str) -> str:
        return OPEN_SEARCH if self._is_open(self._search_values(mode)) else CLOSED_SEARCH

    def _instrument_choice(self, mode: str, current: str | None) -> None:
        key = "analysis-instrument"  # outside param_prefix: the preset values are the parameters
        st.session_state[key] = current
        st.segmented_control(
            "Instrument", list(INSTRUMENT_LABELS), format_func=lambda o: INSTRUMENT_LABELS[o],
            key=key, on_change=self._apply_instrument, args=(key, current, mode),
            help="Sets the precursor and fragment mass tolerances of the search for your mass "
                 "spectrometer. **Orbitrap (HCD)**: 10 ppm and 20 ppm. **timsTOF**: 20 ppm and "
                 "40 ppm. **Ion trap (CID)**: 20 ppm and 0.5 Da.",
        )

    def _apply_instrument(self, key, current, mode):
        chosen = st.session_state.get(key)
        if chosen and chosen != current:
            self._use_instrument(chosen, mode)

    def _use_instrument(self, chosen, mode):
        was_open = self.current_search(mode) == OPEN_SEARCH
        self.ui._apply_preset_callback(chosen)
        if was_open:
            # The instrument presets also return to a standard search
            self.parameter_manager.apply_preset(OPEN_SEARCH)

    def _search_choice(self, mode: str, current_search: str, current_instrument: str | None) -> None:
        key = "analysis-search"
        st.session_state[key] = current_search
        st.segmented_control(
            "Search", list(SEARCH_LABELS), format_func=lambda o: SEARCH_LABELS[o], key=key,
            on_change=self._apply_search, args=(key, current_search, current_instrument),
            help=SEARCH_HELP,
        )

    def _apply_search(self, key, current, current_instrument):
        chosen = st.session_state.get(key)
        if not chosen or chosen == current:
            return
        if chosen == CLOSED_SEARCH and current_instrument:
            # Back to the instrument's own precursor tolerance
            self.ui._apply_preset_callback(current_instrument)
        else:
            self.ui._apply_preset_callback(chosen)

    def _selected_runs(self) -> list:
        selected = self.current_params().get("mzML-files") or []
        return [selected] if isinstance(selected, str) else list(selected)

    def _recommend_instrument(self, current: str | None) -> None:
        """Recommend the instrument preset the selected runs were measured with,
        and offer to apply it when the current settings differ."""
        selected = self._selected_runs()
        if not selected:
            return
        converted_dir = Path(self.workflow_dir, "converted-mzML")
        mtimes = tuple(
            Path(p).stat().st_mtime if Path(p).exists() else None for p in selected
        ) + ((converted_dir.stat().st_mtime,) if converted_dir.exists() else ())
        recommendation = _instrument_recommendation(tuple(selected), str(converted_dir), mtimes)
        if recommendation.preset is None:
            if recommendation.runs:
                st.warning(recommendation.reason(), icon="⚠️")
            return
        text = recommendation.reason()
        if recommendation.unreadable:
            text += " Thermo .raw files are read once they are converted, when the analysis starts."
        if current == recommendation.preset:
            st.success(text, icon="✅")
            return
        c1, c2 = st.columns([4, 1], vertical_alignment="center")
        c1.info(text, icon="💡")
        c2.button(
            f"Use {INSTRUMENT_LABELS[recommendation.preset]} settings",
            key="apply_recommended_instrument",
            on_click=self._use_instrument,
            args=(recommendation.preset, self.mode()),
            use_container_width=True,
        )

    def analysis_section(self) -> None:
        """Your analysis: the runs and the database, the experiment type, the
        instrument and the search, and what a start will do. Drawn above the
        parameters, not in a fragment: the tabs of configure() follow these choices."""
        with st.container(border=True):
            st.markdown("##### Your analysis")
            # The examples are in place before the first selection takes every file
            self.ui._ensure_fallback("mzML-files", EXAMPLE_MZML)
            self.preselect_all_inputs("mzML-files")
            # reactive: drawn in the page itself, so a change reruns the page.
            # Browse opens the runs by folder, with the upload and imports in the same dialog
            self.ui.select_input_file(
                "mzML-files", name="MS data", multiple=True, reactive=True,
                upload=dict(file_types=MS_DATA_TYPES, fallback=EXAMPLE_MZML, suggest_folders=True),
            )
            self.ui.select_input_file(
                "fasta-file", name="Protein database (FASTA)", multiple=False, reactive=True,
                upload=dict(file_types=DATABASE_TYPES, fallback=EXAMPLE_FASTA, uniprot=True),
            )
            self.params = self.parameter_manager.get_parameters_from_json()
            self.apply_decoy_default()

            mode = self._choice(
                "analysis-mode", list(MODE_LABELS), MODE_LABELS, "LFQ", "Experiment type",
                help=" ".join(f"**{MODE_LABELS[m]}**: {h}" for m, h in MODE_HELP.items()),
            )
            current_instrument = self.current_instrument(mode)
            current_search = self.current_search(mode)
            c1, c2 = st.columns([3, 2])
            with c1:
                self._instrument_choice(mode, current_instrument)
            with c2:
                self._search_choice(mode, current_search, current_instrument)
            if current_instrument is None:
                st.caption(
                    "The mass tolerances fit none of the instruments: they are your own "
                    "(Database search tab), or pick your instrument to set them."
                )
            self._recommend_instrument(current_instrument)

            with st.expander("🧬 Reanalyse a PRIDE dataset, set up from its sample metadata (SDRF)"):
                pride_import_section(self, {"mzML-files": list(EXAMPLE_MZML), "fasta-file": EXAMPLE_FASTA})

            params = self.current_params()
            runs = len(self._selected_runs())
            database = params.get("fasta-file")
            st.info(self.summary(mode, runs, current_instrument, current_search,
                                 Path(database).name if isinstance(database, str) else None), icon="🧭")

    def run_bar_summary(self) -> str:
        """One line for the Run bar: analysis, runs and database."""
        params = self.current_params()
        mode = params.get("analysis-mode", "LFQ")
        runs = len(self._selected_runs())
        database = params.get("fasta-file")
        parts = ["Label-free" if mode == "LFQ" else "Isobaric (TMT/iTRAQ)",
                 f"{runs} run{'s' if runs != 1 else ''}"]
        if isinstance(database, str) and database:
            parts.append(Path(database).name)
        return " · ".join(parts) + " · settings are saved when the analysis starts"

    @staticmethod
    def summary(mode: str, runs: int, instrument_preset: str | None, search: str,
                database: str | None) -> str:
        """What a start does, in one paragraph."""
        if not runs:
            return "Choose the runs to analyse with **📂 Browse**."
        what = "label-free" if mode == "LFQ" else "isobaric (TMT/iTRAQ)"
        article = "An" if what[0] in "aeiou" else "A"
        text = (f"**{article} {what} analysis of {runs} run{'s' if runs != 1 else ''}**: "
                f"{'an open' if search == OPEN_SEARCH else 'a'} database search")
        if database:
            text += f" against **{database}** with decoys and common contaminants"
        else:
            text += ", once you choose a protein database"
        text += (f" ({INSTRUMENT_LABELS[instrument_preset]} tolerances)" if instrument_preset
                 else " (your own tolerances)")
        text += ", rescoring with Percolator, a 1 % FDR filter by default, protein inference and "
        text += ("quantification across the runs, then the comparison of your groups."
                 if mode == "LFQ" else
                 "reporter ion quantification, then the comparison of your groups.")
        return text

    @st.fragment
    def configure(self) -> None:
        # The runs, the database and the experiment type are chosen in Your analysis
        self.params = self.parameter_manager.get_parameters_from_json()
        if self.mode() == "LFQ":
            self.render_lfq_tabs()
        else:
            self.render_tmt_tabs()

    def preselect_all_inputs(self, key: str) -> None:
        """Select every uploaded file of ``key`` until the user makes a choice.

        A fresh workspace then runs as it is (all uploaded runs, or the demo
        data); once the selection is saved, even as empty, it is left alone.
        """
        if key in self.params:
            return
        path = Path(self.workflow_dir, "input-files", key)
        if not path.is_dir():
            return
        options = sorted(str(f) for f in path.iterdir() if "external_files.txt" not in str(f))
        if not options:
            return
        params = self.parameter_manager.get_parameters_from_json()
        params[key] = options
        with open(self.parameter_manager.params_file, "w", encoding="utf-8") as f:
            json.dump(params, f, indent=4)
        self.params = params
        self.ui.params = params

    def apply_decoy_default(self) -> None:
        """Default decoy generation to on, unless the selected FASTA already has decoys.

        Runs once per database: a user's choice is kept until another FASTA is
        selected. When decoys are found, their prefix is kept as "decoy-prefix" so
        Percolator and IDFilter recognise them (ProSE detects them itself).
        """
        prefix = self.parameter_manager.param_prefix
        fasta_key = st.session_state.get(f"{prefix}fasta-file", self.params.get("fasta-file"))
        if not fasta_key:
            return
        try:
            fasta = Path(self.file_manager.get_files([fasta_key])[0])
        except (ValueError, IndexError):
            return
        if not fasta.is_file():
            return

        decoy_prefix = detect_decoy_prefix(str(fasta), fasta.stat().st_mtime)
        marker = f"{fasta.name}:{decoy_prefix or ''}"
        if self.params.get("decoy-check") == marker:
            return

        generate = decoy_prefix is None
        params = self.parameter_manager.get_parameters_from_json()
        params["decoy-check"] = marker
        params["generate-decoys"] = generate
        if decoy_prefix:
            params["decoy-prefix"] = decoy_prefix
        with open(self.parameter_manager.params_file, "w", encoding="utf-8") as f:
            json.dump(params, f, indent=4)
        st.session_state[f"{prefix}decoy-check"] = marker
        st.session_state[f"{prefix}generate-decoys"] = generate
        self.params = params
        self.ui.params = params

    @st.fragment
    def render_psm_fdr_widget(self) -> None:
        """PSM FDR as quick choices (0.1 / 1 / 5 %) plus a number field, read back as a q-value.

        The level is stored as "psm-fdr-percent"; the chips only write that field, so a
        value typed in the field and one picked from a chip are the same parameter.
        """
        number_key = f"{self.parameter_manager.param_prefix}psm-fdr-percent"
        chip_key = "psm-fdr-quick-choice"  # outside param_prefix: not a parameter

        def chip_changed():
            choice = st.session_state.get(chip_key)
            if choice is not None:
                st.session_state[number_key] = float(choice)

        def number_changed():
            st.session_state[chip_key] = matching_quick_choice(st.session_state.get(number_key))

        current = st.session_state.get(number_key, self.params.get("psm-fdr-percent", 1.0))
        if chip_key not in st.session_state:
            st.session_state[chip_key] = matching_quick_choice(current)

        st.markdown("**False discovery rate for peptide-spectrum matches (PSMs)**")
        c1, c2, c3 = st.columns([2, 2, 3], vertical_alignment="bottom")
        with c1:
            st.segmented_control(
                "Common choices",
                options=list(FDR_QUICK_CHOICES),
                format_func=percent_label,
                key=chip_key,
                on_change=chip_changed,
                help="0.1 % for a strict list, 1 % is the community standard, "
                     "5 % for exploration.",
            )
        with c2:
            self.ui.input_widget(
                key="psm-fdr-percent",
                default=1.0,
                name="PSM FDR (%)",
                widget_type="number",
                min_value=0.001,
                max_value=100.0,
                step_size=1.0,
                help="Keep PSMs with a q-value up to this level: 1 % keeps q ≤ 0.01. "
                     "100 % switches the filter off and passes every PSM on.",
                on_change=number_changed,
                reactive=True,
            )
        level = st.session_state.get(number_key, current)
        with c3:
            if fdr_applied(level):
                st.markdown(
                    f"<div style='font-size:1.15rem;padding-bottom:0.45rem'>"
                    f"<b>{percent_label(level)}</b> &nbsp;=&nbsp; keep PSMs with "
                    f"<b>{q_label(level)}</b></div>",
                    unsafe_allow_html=True,
                )
            else:
                st.markdown(
                    "<div style='padding-bottom:0.45rem'><b>100 %</b>: PSM FDR filter "
                    "<b>not applied</b>, every PSM is passed on</div>",
                    unsafe_allow_html=True,
                )

    def psm_fdr(self) -> float:
        """PSM FDR threshold as a fraction (1.0 means no filtering)."""
        return min(float(self.params.get("psm-fdr-percent", 1.0)), 100.0) / 100.0

    def filter_psms(self, in_files: list, out_files: list, tool_instance_name: str = "IDFilter") -> bool:
        """Run IDFilter at the configured PSM FDR, or pass PSMs through at 100%."""
        fdr = self.psm_fdr()
        if fdr >= 1.0:
            self.logger.log("FDR level is 100%: skipping PSM FDR filtering")
            for src, dst in zip(in_files, out_files):
                shutil.copy(src, dst)
            return True
        return self.executor.run_topp(
            "IDFilter",
            {"in": in_files, "out": out_files},
            {"score:type_peptide": "q-value", "score:psm": fdr},
            tool_instance_name=tool_instance_name,
        )

    def render_contaminant_option(self):
        """Checkbox for appending the cRAP contaminants to the search database."""
        self.ui.input_widget(
            key="add-contaminants",
            default=True,
            name="Add common contaminants (cRAP)",
            widget_type="checkbox",
            help=(
                "Search the 116 proteins of the common Repository of Adventitious Proteins "
                "(keratins, trypsin, BSA, ...) with your database. Their accessions start with "
                f"{CONTAMINANT_PREFIX}; they are flagged in the PSM tables and left out of the "
                "protein tables and the downstream analysis unless you show them."
            ),
        )

    def render_search_settings(self, instance: str, defaults: dict) -> None:
        """ProSE database search settings (label-free: "ProSE", TMT: "ProSE-TMT")."""
        info = SEARCH_INFO
        if not self.params.get("generate-decoys", True):
            info += SEARCH_DECOY_INFO.format(prefix=self.params.get("decoy-prefix", "rev_"))
        st.info(info)
        self.ui.input_TOPP(
            "ProSE",
            custom_defaults=defaults,
            include_parameters=SEARCH_INCLUDE,
            exclude_parameters=SEARCH_EXCLUDE,
            tool_instance_name=instance,
        )

    def search_tolerance(self, instance: str) -> tuple:
        """Fragment tolerance of the search and whether it is in ppm, for the spectrum views."""
        params = self.parameter_manager.get_topp_parameters("ProSE", tool_instance_name=instance)
        tolerance = float(params.get("Search:fragment:mass_tolerance", 20.0))
        return tolerance, params.get("Search:fragment:mass_tolerance_unit", "ppm") == "ppm"

    def run_search(self, instance: str, in_mzml: list, out_idxml: list, database_fasta: Path) -> bool:
        """Search the spectra with ProSE against the database (with its decoys).

        ProSE searches several files in one process: it builds the fragment
        index once and reuses it, and each file gets all the threads of the
        process. One process per file built the index once per file and,
        with many files, ran as many processes as there are threads, each
        with one thread and an index of its own in memory. The files are
        split into ``search_batches()``: one process per
        SEARCH_THREADS_PER_PROCESS threads of the budget.
        """
        batches = search_batches(len(in_mzml), self.executor._get_max_threads())
        self.logger.log(
            f"Searching {len(in_mzml)} file(s) in {len(batches)} ProSE process(es)", 1
        )
        return self.executor.run_topp(
            "ProSE",
            {
                "in": [[in_mzml[i] for i in b] for b in batches],
                "out_idxml": [[out_idxml[i] for i in b] for b in batches],
            },
            {"database": str(database_fasta), "Search:decoys": "auto"},
            tool_instance_name=instance,
        )

    def render_lfq_tabs(self):
        st.subheader("Label-free experiment")
        t = st.tabs([
            "**1 · Samples and groups**", "**2 · Database search**", "**3 · Rescoring**",
            "**4 · False discovery rate**", "**5 · Quantification**", "**6 · Spectral library (optional)**",
        ])

        with t[0]:
            self.render_lfq_groups()

        with t[1]:
            # Checkbox for decoy generation
            # reactive=True ensures the parent configure() fragment re-runs when checkbox changes,
            # so conditional UI (DecoyDatabase settings) updates immediately
            self.ui.input_widget(
                key="generate-decoys",
                default=True,
                name="Generate Decoy Database",
                widget_type="checkbox",
                help="Generate decoy sequences for FDR calculation. Switched off automatically when the selected FASTA already contains decoys (e.g. DECOY_ or rev_ accessions).",
                reactive=True,
            )
            self.render_contaminant_option()

            # Reload params to get current checkbox value after it was saved
            self.params = self.parameter_manager.get_parameters_from_json()

            # Show DecoyDatabase settings if generating decoys
            if self.params.get("generate-decoys", True):
                st.info(DECOY_INFO)
                self.ui.input_TOPP(
                    "DecoyDatabase",
                    custom_defaults={
                        "decoy_string": "rev_",
                        "decoy_string_position": "prefix",
                        "method": "reverse",
                    },
                    include_parameters=["method"],
                )

            self.render_search_settings("ProSE", SEARCH_DEFAULTS)

        with t[2]:
            st.info(PERCOLATOR_INFO)
            # decoy_pattern is always derived from upstream, never shown
            percolator_include = ["post_processing_tdc", "score_type", "subset_max_train"]

            self.ui.input_TOPP(
                "PercolatorAdapter",
                custom_defaults={
                    "threads": 8,
                    "subset_max_train": 300000,
                    "decoy_pattern": "rev_",
                    "score_type": "pep",
                    "post_processing_tdc": True,
                },
                include_parameters=percolator_include,
                flag_parameters=["post_processing_tdc"],
                exclude_parameters=["out_type"],
            )
            self.render_prediction_features_widget()

        with t[3]:
            st.info(FDR_INFO)
            self.render_psm_fdr_widget()
            self.ui.input_TOPP(
                "IDFilter",
                custom_defaults={
                    "threads": 2,
                },
                exclude_parameters=["type_protein", "type_peptide", "score:psm"],
            )

        with t[5]:
            st.info("""
            **Optional.** Builds a spectral library from the identified peptides, for analysing
            data-independent acquisition (DIA/SWATH) runs of the same samples later, e.g. in
            OpenSWATH. Not needed for the protein comparison; leave it off if you do not plan a
            DIA experiment.
            """)

            self.ui.input_widget(
                key="generate-library",
                default=False,
                name="Generate Spectral Library",
                widget_type="checkbox",
                help="Enable spectral library generation using EasyPQP.",
                reactive=True,
            )
            self.params = self.parameter_manager.get_parameters_from_json()

            if self.params.get("generate-library", False):
                # FDR options
                self.ui.input_widget(
                    key="library-use-fdr",
                    default=False,
                    name="Apply additional FDR filtering",
                    widget_type="checkbox",
                    help="If disabled (recommended), uses --nofdr since IDFilter already applied FDR.",
                    reactive=True,
                )
                self.params = self.parameter_manager.get_parameters_from_json()

                if self.params.get("library-use-fdr", False):
                    self.ui.input_widget(
                        key="library-psm-fdr",
                        default=0.01,
                        name="PSM FDR Threshold",
                        widget_type="number",
                        min_value=0.001,
                        max_value=1.0,
                        step_size=0.01,
                    )

                # Decoy options
                self.ui.input_widget(
                    key="library-generate-decoys",
                    default=True,
                    name="Generate Decoy Library",
                    widget_type="checkbox",
                    reactive=True,
                )
                self.params = self.parameter_manager.get_parameters_from_json()

                if self.params.get("library-generate-decoys", True):
                    self.ui.input_widget(
                        key="library-decoy-method",
                        default="shuffle",
                        name="Decoy Method",
                        widget_type="selectbox",
                        options=["shuffle", "reverse", "pseudo-reverse", "shift"],
                    )

        with t[4]:
            st.info("""
            **How proteins are quantified.** ProteomicsLFQ measures each identified peptide's
            signal over its elution time in every run, aligns the runs so the same peptide is
            matched across them, and adds up the peptides of each protein. Only peptides unique
            to one protein are used, so shared peptides cannot blur the comparison.

            * **intThreshold**: the faintest signal that counts as a peptide (raise it for noisy data).
            * **psmFDR** / **proteinFDR**: the false discovery rate at peptide and protein level
              (0.01 = 1 %, the community standard).
            * **Match between runs**: a peptide identified in one run is also quantified in the
              runs where its signal was measured but not identified. Fewer missing values, at the
              cost of a longer run and some transfers to the wrong signal.
            """)
            self.ui.input_widget(
                key=MATCH_BETWEEN_RUNS_KEY,
                default=False,
                name="Match between runs",
                widget_type="checkbox",
                help=(
                    "Also quantify a peptide in runs where it was not identified, using the "
                    "aligned retention time and m/z of its identification in another run. "
                    "Fills in missing values; best for runs of similar samples on one LC setup."
                ),
            )
            self.ui.input_TOPP(
                "ProteomicsLFQ",
                custom_defaults={
                    "threads": 12,
                    "feature_with_id_min_score": 0.1,
                    "Seeding:intThreshold": 1000.0,
                    "psmFDR": 0.01,
                    "proteinFDR": 0.01,
                    "picked_proteinFDR": "true",
                    "alignment_order": "star",
                    "protein_quantification": "unique_peptides",
                    "quantification_method": "feature_intensity",
                    "protein_inference": "aggregation",
                },
                include_parameters=["intThreshold", "psmFDR", "proteinFDR"],
            )

    def render_lfq_groups(self):
        """Experimental design: group, biological replicate and fraction of every file."""
        st.markdown(
            "Give every file the name of its experimental group (**Condition**), for example "
            "**control** and **treated**. Files of one condition are its replicates, and the group "
            "comparison tests every protein between the conditions. Use at least two conditions, "
            "ideally with three or more biological replicates each."
        )
        st.markdown(
            "If a sample was **fractionated** before the MS run (for example by high-pH reversed "
            "phase), give its files the same condition and replicate and number them as fractions "
            "1, 2, 3, ...; their signals are then combined into one sample. Unfractionated samples "
            "keep fraction 1. You can also import a design.tsv or an SDRF file."
        )

        mzml_keys = self.params.get("mzML-files")

        if not mzml_keys:
            st.warning("No mzML files available. Please upload mzML files first.")
            return

        try:
            mzml_files = self.file_manager.get_files(mzml_keys)
        except ValueError:
            st.warning("Selected mzML files are not available.")
            return

        # Current mzML filenames
        current_filenames = {Path(mz).name for mz in mzml_files}

        # The table edits the design file; the groups stay in the mzML-group-<file>
        # parameters, which the results pages read. Bring the table up to date with
        # them first, as they can also be set elsewhere (PRIDE import).
        files = st.session_state.get(
            f"{self.parameter_manager.param_prefix}mzML-files", mzml_keys
        )
        if isinstance(files, str):
            files = [files]
        path = design_path(self.workflow_dir)
        saved = path.read_text(encoding="utf-8") if path.is_file() else None
        stored = load_design(self.workflow_dir, files)
        groups = initial_groups(stored["file"], self.params, demo_group)
        seeded = seed_design(stored, groups, renumber=saved is None)
        if to_openms_design(seeded, files) != saved:
            save_design(self.workflow_dir, seeded, files)
            # A new editor key, so the table shows the file instead of earlier edits
            version_key = f"_design-{Path(self.workflow_dir).name}-mzML-files"
            st.session_state[version_key] = st.session_state.get(version_key, 0) + 1

        self.ui.experimental_design("mzML-files", name="Samples")

        prefix = self.parameter_manager.param_prefix
        changed = False
        for filename, group in groups_from_design(load_design(self.workflow_dir, files)).items():
            if self.params.get(f"mzML-group-{filename}") != group:
                st.session_state[f"{prefix}mzML-group-{filename}"] = group
                changed = True
        if changed:
            self.parameter_manager.save_parameters()

        # Reload params to get current values
        self.params = self.parameter_manager.get_parameters_from_json()

        # Clean up orphaned group params from previously selected files
        orphaned_keys = [
            k for k in self.params.keys()
            if k.startswith("mzML-group-") and k[11:] not in current_filenames
        ]
        for key in orphaned_keys:
            del self.params[key]
        if orphaned_keys:
            self.parameter_manager.save_parameters()

    def render_tmt_tabs(self):
        st.subheader("TMT / iTRAQ experiment")
        # Create tabs for different analysis steps.
        t = st.tabs([
            "**1 · Samples and groups**", "**2 · Labels and decoys**", "**3 · Database search**",
            "**4 · Rescoring**", "**5 · False discovery rate**", "**6 · Proteins and quantification**",
        ])

        with t[0]:
            self.render_tmt_groups()
        with t[1]:
            # Checkbox for decoy generation
            # reactive=True ensures the parent configure() fragment re-runs when checkbox changes,
            # so conditional UI (DecoyDatabase settings) updates immediately
            self.ui.input_widget(
                key="generate-decoys",
                default=True,
                name="Generate Decoy Database",
                widget_type="checkbox",
                help="Generate decoy sequences for FDR calculation. Switched off automatically when the selected FASTA already contains decoys (e.g. DECOY_ or rev_ accessions).",
                reactive=True,
            )
            self.render_contaminant_option()

            # Reload params to get current checkbox value after it was saved
            self.params = self.parameter_manager.get_parameters_from_json()

            # Show DecoyDatabase settings if generating decoys
            if self.params.get("generate-decoys", True):
                st.info(DECOY_INFO)
                self.ui.input_TOPP(
                    "DecoyDatabase",
                    custom_defaults={
                        "decoy_string": "rev_",
                        "decoy_string_position": "prefix",
                        "method": "reverse",
                    },
                    include_parameters=["method"],
                )

            st.info("""
            **Label type.** Choose the TMT or iTRAQ kit you used (e.g. *tmt11plex*, *tmt16plex*).
            IsobaricAnalyzer reads each channel's reporter-ion intensity from the fragment spectra;
            the defaults fit current Orbitrap TMT methods.
            """)
            self.ui.input_TOPP(
                "IsobaricAnalyzer",
                custom_defaults={
                    "tmt11plex:reference_channel": 126,
                    "type": "tmt11plex",
                    "extraction:select_activation": "auto",
                    "extraction:reporter_mass_shift": 0.002,
                    "extraction:min_reporter_intensity": 0.0,
                    "extraction:min_precursor_purity": 0.0,
                    "extraction:precursor_isotope_deviation": 10.0,
                    "quantification:isotope_correction": "false",
                },
                tool_instance_name="IsobaricAnalyzer-TMT",
                reactive=True,
            )
        with t[2]:
            self.render_search_settings("ProSE-TMT", SEARCH_DEFAULTS_TMT)
        with t[3]:
            st.info(PERCOLATOR_INFO)
            self.ui.input_TOPP(
                "PercolatorAdapter",
                custom_defaults={
                    "subset_max_train": 300000,
                    "decoy_pattern": "DECOY_",
                    "score_type": "pep",
                    "post_processing_tdc": True,
                    "debug": 0,
                },
                flag_parameters=["post_processing_tdc"],
                tool_instance_name="PercolatorAdapter-TMT",
            )
            self.render_prediction_features_widget()

        with t[4]:
            st.info(FDR_INFO)
            self.render_psm_fdr_widget()
            self.ui.input_TOPP(
                "IDFilter",
                exclude_parameters=["type_peptide", "score:psm"],
                tool_instance_name="IDFilter-strict",
            )
        with t[5]:
            st.info(
                "These steps turn the identified peptides into one abundance per protein and channel. "
                "The defaults follow the quantms TMT workflow; most experiments never change them."
            )
            st.markdown("**Map identifications to the reporter intensities** (IDMapper)")
            self.ui.input_TOPP(
                "IDMapper",
                custom_defaults={
                    "threads": 8,
                    "debug": 0,
                },
                tool_instance_name="IDMapper-TMT",
            )
        with t[5]:
            st.markdown("**Combine the runs** (FileMerger)")
            self.ui.input_TOPP(
                "FileMerger",
                custom_defaults={
                    "in_type": "consensusXML",
                    "append_method": "append_cols",
                    "annotate_file_origin": True,
                    "threads": 8,
                },
                flag_parameters=["annotate_file_origin"],
                tool_instance_name="FileMerger-TMT",
            )
        with t[5]:
            st.markdown("**Infer proteins from peptides** (ProteinInference)")
            self.ui.input_TOPP(
                "ProteinInference",
                custom_defaults={
                    "threads": 8,
                    "picked_decoy_string": "DECOY_",
                    "picked_fdr": "true",
                    "protein_fdr": "true",
                    "Algorithm:use_shared_peptides": "true",
                    "Algorithm:annotate_indistinguishable_groups": "true",
                    "Algorithm:score_type": "PEP",
                    "Algorithm:score_aggregation_method": "best",
                    "Algorithm:min_peptides_per_protein": 1,
                },
                tool_instance_name="ProteinInference-TMT",
            )
        with t[5]:
            st.markdown("**Protein-level false discovery rate** (IDFilter)")
            # A single checkbox widget for workflow logic.
            # self.ui.input_widget("run-python-script", False, "Run custom Python script") *
            # Generate input widgets for a custom Python tool, located at src/python-tools.
            # Parameters are specified within the file in the DEFAULTS dictionary.
            # self.ui.input_python("example") *
            self.ui.input_TOPP(
                "IDFilter",
                custom_defaults={
                    "score:type_protein": "q-value",
                    "score:proteingroup": 0.01,
                    "score:psm": 0.01,
                    "delete_unreferenced_peptide_hits": True,
                    "remove_decoys": True
                },
                flag_parameters=["delete_unreferenced_peptide_hits", "remove_decoys"],
                tool_instance_name="IDFilter-lenient",
            )
        with t[5]:
            st.markdown("**Resolve peptides matched to several features** (IDConflictResolver)")
            self.ui.input_TOPP(
                "IDConflictResolver",
                custom_defaults={
                    "threads": 4,
                },
                tool_instance_name="IDConflictResolver-TMT",
            )

        with t[5]:
            st.markdown("**Sum peptides into protein abundances** (ProteinQuantifier)")
            self.ui.input_TOPP(
                "ProteinQuantifier",
                custom_defaults={
                    "method": "top",
                    "top:N": 3,
                    "top:aggregate": "median",
                    "top:include_all": True,
                    "threads": 8,
                    "debug": 0,
                },
                flag_parameters=["top:include_all"],
                tool_instance_name="ProteinQuantifier-TMT",
            )
    def render_tmt_groups(self):
        """Per-channel group names for the label type chosen in IsobaricAnalyzer."""
        st.markdown(
            "Give every label channel the name of the experimental group of the sample in it, for "
            "example **control** and **treated**. Channels with the same name are replicates. "
            "Type **skip** for an empty channel or a reference/pool channel you do not want compared. "
            "The number of channels follows the label type chosen in tab 2."
        )

        latest_params = self.parameter_manager.get_parameters_from_json()
        type_key = (
            f"{self.parameter_manager.topp_param_prefix}"
            "IsobaricAnalyzer-TMT:1:type"
        )
        selected_type = str(
            st.session_state.get(type_key)
            or latest_params.get("IsobaricAnalyzer-TMT", {}).get("type")
            or "tmt11plex"
        ).lower()

        m = re.search(r'\d+', selected_type)
        is_supported_type = any(label in selected_type for label in ["tmt", "itraq"])
        if not m or not is_supported_type:
            st.warning("Choose a TMT or iTRAQ label type in tab 2 (Labels and decoys) first.")
        else:
            num_plex = int(m.group())
            channels = [f"sample{i+1}" for i in range(num_plex)]
            reporter = ISOBARIC_CHANNELS.get(selected_type, [])

            def tags(idx):
                return f" ({reporter[idx]})" if idx < len(reporter) else ""
            st.caption(f"Isobaric type: **{selected_type}** - {num_plex} channels")

            for row_start in range(0, num_plex, 2):
                c1, c2 = st.columns(2)

                left_idx = row_start
                left_channel = channels[left_idx]
                with c1:
                    self.ui.input_widget(
                        key=f"TMT-group-{left_channel}",
                        default="",
                        name=f"Group for channel {left_idx + 1}{tags(left_idx)}",
                        widget_type="text",
                        help="e.g. control, case, skip",
                    )

                right_idx = row_start + 1
                if right_idx < num_plex:
                    right_channel = channels[right_idx]
                    with c2:
                        self.ui.input_widget(
                            key=f"TMT-group-{right_channel}",
                            default="",
                            name=f"Group for channel {right_idx + 1}{tags(right_idx)}",
                            widget_type="text",
                            help="e.g. control, case, skip",
                        )

            # Remove orphaned params from a previously selected larger plex
            self.params = self.parameter_manager.get_parameters_from_json()
            valid_keys = {f"TMT-group-{ch}" for ch in channels}
            orphaned = [k for k in self.params if k.startswith("TMT-group-") and k not in valid_keys]
            if orphaned:
                for k in orphaned:
                    del self.params[k]
                self.parameter_manager.save_parameters()

    def render_prediction_features_widget(self) -> None:
        """Opt-in switch: AlphaPeptDeep features for Percolator (shared by both modes)."""
        st.markdown("**Rescore with predicted spectra and retention times** (AlphaPeptDeep)")
        available = peptdeep_available()
        self.ui.input_widget(
            key="rescore-predictions",
            default=False,
            name="Add predicted-spectrum features",
            widget_type="checkbox",
            help=(
                "Predicts every match's fragment intensities and retention time with AlphaPeptDeep "
                "and gives Percolator three more features: spectral angle and cosine similarity "
                "between the measured and the predicted spectrum, and the retention-time "
                "difference. Correct matches look like their prediction; wrong ones rarely do, "
                "so more matches pass the same FDR. On a CPU it adds about 15 s to load the models "
                "plus roughly 4 s per 1,000 matches."
            ),
            reactive=True,
        )
        if not available:
            st.caption("AlphaPeptDeep (the peptdeep package) is not installed here, so this option has no effect.")
            return
        self.params = self.parameter_manager.get_parameters_from_json()
        if self.params.get("rescore-predictions", False):
            self.ui.input_widget(
                key="rescore-instrument",
                default="auto",
                name="Instrument model",
                widget_type="selectbox",
                options=["auto", "QE", "Lumos", "Exploris", "Astral", "timsTOF", "SciexTOF", "Velos", "Elite"],
                help=(
                    "The instrument group whose fragmentation AlphaPeptDeep predicts. **auto** reads "
                    "it from each mzML header and falls back to Lumos. The collision energy is "
                    "calibrated on each file's confident matches."
                ),
            )

    def add_prediction_features(self, search_results: list, in_mzML: list, search_instance: str):
        """Annotate the search results with AlphaPeptDeep features when the option is on.

        Returns the idXMLs Percolator should read: the search results unchanged
        when the option is off, the annotated copies when it is on, and None when
        it is on but failed (the run stops).
        """
        if not self.params.get("rescore-predictions", False):
            return search_results
        if not peptdeep_available():
            self.logger.log("❌ Predicted-spectrum features need the peptdeep package, which is not installed.")
            self.logger.log("Workflow stopped due to error")
            return None
        out_dir = Path(self.results_dir, "prediction_features")
        out_dir.mkdir(parents=True, exist_ok=True)
        annotated = [str(out_dir / f"{Path(f).stem}_pd.idXML") for f in search_results]
        fragment_tolerance, fragment_ppm = self.search_tolerance(search_instance)
        params_file = out_dir / "params.json"
        summary_file = out_dir / "summary.json"
        params_file.write_text(json.dumps({
            "in": search_results,
            "mzml": in_mzML,
            "out": annotated,
            "summary": str(summary_file),
            "fragment_tolerance": fragment_tolerance,
            "fragment_ppm": fragment_ppm,
            "instrument": self.params.get("rescore-instrument", "auto"),
        }, indent=2), encoding="utf-8")
        self.logger.log("🧠 Predicting spectra and retention times (AlphaPeptDeep)...")
        script = Path("src", "python-tools", "peptdeep_rescore.py")
        ok = self.executor.run_command([sys.executable, str(script), str(params_file)])
        params_file.unlink(missing_ok=True)
        if not ok or not all(Path(f).exists() for f in annotated):
            self.logger.log("Workflow stopped due to error")
            return None
        for s in json.loads(summary_file.read_text(encoding="utf-8")):
            self.logger.log(
                f"  {s['file']}: {s['instrument']}, NCE {s['nce']:g}, {s['predicted']} of {s['psms']} matches "
                f"predicted; median spectral angle targets {s['median_angle_targets'] or 0:.2f}, "
                f"decoys {s['median_angle_decoys'] or 0:.2f}"
            )
        self.logger.log("✅ Prediction features added")
        return annotated

    def _has_identifications(self, idxml_files: list, stage: str) -> bool:
        """Return False, and log why, when no idXML file holds a peptide hit.

        Every later step (rescoring, filtering, quantification and the result
        views) needs at least one identification, so an empty search or an
        FDR filter that removes everything has to stop the run with a message
        the user can act on instead of crashing further down.
        """
        n_hits = sum(
            len(pep.getHits())
            for f in idxml_files
            if Path(f).exists()
            for pep in load_idxml(f)[1]
        )
        if n_hits > 0:
            return True
        if stage == "FDR filtering":
            reason = (
                f" at PSM FDR {fdr_label(self.psm_fdr() * 100)}. The threshold may be "
                "stricter than any match reached (try 1 % or 5 %), there may be too few "
                "spectra to estimate an FDR, or most matches are wrong: check that the "
                "FASTA database matches the sample species and that the tolerances and "
                "enzyme fit the data."
            )
        else:
            reason = (
                ". Check that the FASTA database matches the sample species, that the "
                "precursor and fragment tolerances and enzyme fit the data, and that the "
                "files hold enough MS2 spectra."
            )
        self.logger.log(f"❌ No peptide identifications left after {stage}{reason}")
        self.logger.log("Workflow stopped due to error")
        return False

    def write_qpx(self) -> None:
        """Export the filtered PSMs and the run's provenance as QPX Parquet (an export; failures are logged)."""
        self.logger.log("📦 Writing the QPX export (PSMs and provenance)...")
        try:
            log = Path(self.workflow_dir, "logs", "commands-and-run-times.log")
            out = export_run(
                self.results_dir, self.workflow_dir, self.params,
                log.read_text(encoding="utf-8", errors="ignore") if log.is_file() else "",
                app_version(),
            )
            self.logger.log(f"✅ QPX export: {out['psms']} PSMs in {Path(out['psm']).name}")
        except Exception as e:  # an export must not fail the run
            self.logger.log(f"⚠️ QPX export failed: {e}")

    def write_identified_spectra(self, mzml_files: list, idxml_files: list) -> None:
        """Keep the spectra of the final PSMs with the run (``results/spectra``).

        A run opened from its ZIP shows its PSM spectra from them when the
        input mzML files are not there. A failure only costs that, so it is
        logged and the run goes on.
        """
        try:
            written = write_identified_spectra(mzml_files, idxml_files, self.results_dir)
            self.logger.log(
                f"🧾 Kept {sum(written.values())} identified spectra in {SPECTRA_DIR}/ for viewing without the mzML files"
            )
        except Exception as error:  # noqa: BLE001 - the results are complete without them
            shutil.rmtree(Path(self.results_dir, SPECTRA_DIR), ignore_errors=True)
            self.logger.log(f"⚠️ Could not keep the identified spectra with the run: {error}")

    def prepare_search_database(self, fasta_path: Path):
        """Build the database to search: the FASTA, the cRAP contaminants and decoys.

        Returns ``(database_fasta, decoy_string)``, or None when a TOPP tool failed.
        """
        # Contaminants (cRAP with CON__ accessions) are searched with the targets
        add_contaminants = self.params.get("add-contaminants", True)
        db_dir = Path(self.results_dir, "database")
        con_fasta = db_dir / "contaminants.fasta"
        if add_contaminants:
            n_con, n_in_db = write_contaminants(con_fasta, targets=fasta_path)
            note = f"; {n_in_db} already in your database are kept as targets" if n_in_db else ""
            self.logger.log(
                f"🧫 Adding {n_con} common lab contaminants (cRAP: keratins, trypsin, ...; "
                f"{CONTAMINANT_PREFIX} accessions, hidden from the results){note}"
            )

        if self.params.get("generate-decoys", True):
            # Get decoy_string from DecoyDatabase params
            decoy_string = self.params.get("DecoyDatabase", {}).get("decoy_string", "rev_")
            if add_contaminants:
                # DecoyDatabase reads both files, so the contaminants get decoys too
                decoy_fasta = db_dir / f"{fasta_path.stem}.con.decoy.fasta"
                decoy_inputs = [str(fasta_path), str(con_fasta)]
                stale = decoy_fasta.exists() and decoy_fasta.stat().st_mtime < max(
                    fasta_path.stat().st_mtime, con_fasta.stat().st_mtime
                )
            else:
                decoy_fasta = fasta_path.with_suffix(".decoy.fasta")
                decoy_inputs = [str(fasta_path)]
                stale = False

            if not decoy_fasta.exists() or stale:
                self.logger.log("🧬 Generating decoy database...")
                st.info("Generating decoy FASTA database...")
                db_dir.mkdir(parents=True, exist_ok=True)
                if not self.executor.run_topp(
                    "DecoyDatabase",
                    {"in": [decoy_inputs], "out": [str(decoy_fasta)]},
                ):
                    self.logger.log("Workflow stopped due to error")
                    return None
                self.logger.log("✅ Decoy database ready")
            st.success(f"Using decoy FASTA: {decoy_fasta.name}")
            database_fasta = decoy_fasta
        else:
            # Prefix detected in the FASTA (workspaces from before: the search engine's setting)
            legacy = "CometAdapter-TMT" if self.params.get("analysis-mode", "LFQ") == "TMT" else "CometAdapter"
            decoy_string = self.params.get("decoy-prefix") or self.params.get(legacy, {}).get(
                "PeptideIndexing:decoy_string", "rev_"
            )
            self.logger.log("📄 Using existing FASTA database")
            st.info(f"Using original FASTA: {fasta_path.name}")
            database_fasta = fasta_path
            if add_contaminants:
                # The FASTA has its own decoys; give the contaminants decoys in the same style.
                # Own instance name, so stored DecoyDatabase settings don't repeat these options.
                con_decoy_fasta = db_dir / "contaminants.decoy.fasta"
                if not self.executor.run_topp(
                    "DecoyDatabase",
                    {"in": [str(con_fasta)], "out": [str(con_decoy_fasta)]},
                    {"decoy_string": decoy_string, "decoy_string_position": "prefix"},
                    tool_instance_name="DecoyDatabase-contaminants",
                ):
                    self.logger.log("Workflow stopped due to error")
                    return None
                database_fasta = db_dir / f"{fasta_path.stem}.con.fasta"
                append_fasta([fasta_path, con_decoy_fasta], database_fasta)
                st.info(f"Added contaminants: {database_fasta.name}")
        return database_fasta, decoy_string

    def alias_converted_runs(self, inputs: list, mzmls: list) -> None:
        """Give each converted run the group of its vendor file.

        The tools and the results pages name a run by its mzML file, while the
        groups were set for the file the user selected (``run.raw``, ``run.d``).
        Written to this run's parameter file, which the results pages read.
        """
        params = self.parameter_manager.get_parameters_from_json() or dict(self.params)
        changed = False
        for original, mzml in zip(inputs, mzmls):
            source, target = f"mzML-group-{Path(original).name}", f"mzML-group-{Path(mzml).name}"
            if source != target and source in self.params:
                params[target] = self.params[source]
                self.params[target] = self.params[source]
                changed = True
        if changed:
            with open(self.parameter_manager.params_file, "w", encoding="utf-8") as f:
                json.dump(params, f, indent=4)

    def execution(self) -> bool:
        """
        Refactored TOPP workflow execution:
        - Per-sample: ProSE -> PercolatorAdapter -> IDFilter
        - Cross-sample: ProteomicsLFQ (single combined output)
        """
        # ================================
        # 0️⃣ Input validation
        # ================================
        if not self.params.get("mzML-files"):
            st.error("No mzML files selected.")
            return False

        if not self.params.get("fasta-file"):
            st.error("No FASTA file selected.")
            return False

        in_mzML = self.file_manager.get_files(self.params["mzML-files"])
        fasta_file = self.file_manager.get_files([self.params["fasta-file"]])[0]

        # Thermo .raw files and Bruker .d folders are converted to mzML first
        converted = self.executor.convert_to_mzml(in_mzML)
        if converted is None:
            self.logger.log("ERROR: A run could not be converted to mzML; see the log above.", 0)
            return False
        if converted != list(in_mzML):
            self.alias_converted_runs(in_mzML, converted)
        in_mzML = converted

        if len(in_mzML) < 1:
            st.error("At least one mzML file is required.")
            return False
        
        fasta_path = Path(fasta_file)

        self.logger.log(f"📂 Loaded {len(in_mzML)} sample(s)")

        prepared = self.prepare_search_database(fasta_path)
        if prepared is None:
            return False
        database_fasta, decoy_string = prepared

        current_mode = self.params.get("analysis-mode", "LFQ")
        st.write(f"Current analysis mode: **{current_mode}**")

        if current_mode == "LFQ":
            self.logger.log("⚙️ Running LFQ workflow")

            # ================================
            # 1️⃣ Directory setup
            # ================================
            results_dir = Path(self.results_dir)
            search_dir = results_dir / "search_results"
            perc_dir = results_dir / "percolator_results"
            filter_dir = results_dir / "psm_filter"
            quant_dir = results_dir / "quant_results"

            results_dir = Path(self.workflow_dir, "input-files")

            for d in [search_dir, perc_dir, filter_dir, quant_dir]:
                d.mkdir(parents=True, exist_ok=True)

            self.logger.log("📁 Output directories created")

            # ================================
            # 2️⃣ File path definitions (per sample)
            # ================================
            search_results = []
            percolator_results = []
            filter_results = []

            for mz in in_mzML:
                stem = Path(mz).stem
                search_results.append(str(search_dir / f"{stem}_prose.idXML"))
                percolator_results.append(str(perc_dir / f"{stem}_per.idXML"))
                filter_results.append(str(filter_dir / f"{stem}_filter.idXML"))

            # ================================
            # 3️⃣ Per-file processing
            # ================================
            for i, mz in enumerate(in_mzML):
                stem = Path(mz).stem
                st.info(f"Processing sample: {stem}")

            self.logger.log("🔬 Starting per-sample processing...")

            # --- ProSE database search ---
            self.logger.log("🔎 Running peptide search...")
            with st.spinner("ProSE database search"):
                if not self.run_search("ProSE", in_mzML, search_results, database_fasta):
                    self.logger.log("Workflow stopped due to error")
                    return False

            # Fragment tolerance of the search, for the spectrum views
            frag_tol, frag_tol_is_ppm = self.search_tolerance("ProSE")

            if not self._has_identifications(search_results, "peptide search"):
                return False

            # Build visualization cache for the search results
            results_dir_path = Path(self.results_dir)
            cache_dir = results_dir_path / "insight_cache"
            cache_dir.mkdir(parents=True, exist_ok=True)

            # Spectra cache (once, shared by all stages): one parquet per input
            # mzML, read only for the idXML of that file, so memory follows one
            # run instead of the whole study. Every input mzML gets a fixed file
            # index, so each PSM links to its own run's spectra.
            filename_to_index = {Path(m).name: i for i, m in enumerate(in_mzML)}
            spectra_dir = Path(self.workflow_dir, "spectra-cache")
            shutil.rmtree(spectra_dir, ignore_errors=True)
            spectra_df = write_spectra_cache(spectra_dir, filename_to_index, in_mzML)

            for idxml_file in search_results:
                idxml_path = Path(idxml_file)

                # Table, Heatmap, SequenceView and LinePlot caches
                if not build_psm_views(
                    idxml_path, cache_dir, spectra_df, filename_to_index, frag_tol, frag_tol_is_ppm,
                    export_dir=None,
                ):
                    self.logger.log(f"⚠️ No identifications in {idxml_path.name}, skipping its result views")

            self.logger.log("✅ Peptide search complete")

            # --- Optional AlphaPeptDeep features ---
            percolator_in = self.add_prediction_features(search_results, in_mzML, "ProSE")
            if percolator_in is None:
                return False

            # --- PercolatorAdapter ---
            self.logger.log("📊 Running rescoring...")
            with st.spinner(f"PercolatorAdapter ({stem})"):
                if not self.executor.run_topp(
                    "PercolatorAdapter",
                    {
                        "in": percolator_in,
                        "out": percolator_results,
                    },
                    {"decoy_pattern": decoy_string},  # Always propagated from upstream
                ):
                    self.logger.log("Workflow stopped due to error")
                    return False

            if not self._has_identifications(percolator_results, "rescoring"):
                return False

            # Build visualization cache for Percolator results
            for idxml_file in percolator_results:
                idxml_path = Path(idxml_file)

                # Table, Heatmap, SequenceView and LinePlot caches
                if not build_psm_views(
                    idxml_path, cache_dir, spectra_df, filename_to_index, frag_tol, frag_tol_is_ppm,
                    export_dir=None,
                ):
                    self.logger.log(f"⚠️ No identifications in {idxml_path.name}, skipping its result views")

            self.logger.log("✅ Rescoring complete")

            # if not Path(percolator_results[i]).exists():
            #     st.error(f"PercolatorAdapter failed for {stem}")
            #     st.stop()

            # --- IDFilter ---
            self.logger.log("🔧 Filtering identifications...")
            with st.spinner(f"IDFilter ({stem})"):
                if not self.filter_psms(percolator_results, filter_results):
                    self.logger.log("Workflow stopped due to error")
                    return False

            if not self._has_identifications(filter_results, "FDR filtering"):
                return False

            # Build visualization cache for Filter results
            for idxml_file in filter_results:
                idxml_path = Path(idxml_file)

                # Table, Heatmap, SequenceView and LinePlot caches and the TSV exports
                if not build_psm_views(
                    idxml_path, cache_dir, spectra_df, filename_to_index, frag_tol, frag_tol_is_ppm,
                    export_dir=results_dir_path / "exports",
                ):
                    self.logger.log(f"⚠️ No identifications in {idxml_path.name}, skipping its result views")
            shutil.rmtree(spectra_dir, ignore_errors=True)

            self.logger.log("✅ Filtering complete")

            # if not Path(filter_results[i]).exists():
            #     st.error(f"IDFilter failed for {stem}")
            #     st.stop()

            # ================================
            # EasyPQP Spectral Library Generation (optional)
            # ================================
            if self.params.get("generate-library", False):
                self.logger.log("📚 Building spectral library with EasyPQP...")
                st.info("Building spectral library with EasyPQP...")
                library_dir = Path(self.results_dir, "library")
                library_dir.mkdir(parents=True, exist_ok=True)

                psms_files, peaks_files = [], []

                for filter_idxml in filter_results:
                    original_stem = Path(filter_idxml).stem.replace("_filter", "")
                    matching_mzml = next((m for m in in_mzML if Path(m).stem == original_stem), None)
                    if not matching_mzml:
                        self.logger.log(f"Warning: No matching mzML found for {filter_idxml}")
                        continue

                    # easypqp library requires specific extensions for file recognition:
                    # - PSM files must contain 'psmpkl' → use .psmpkl extension
                    # - Peak files must contain 'peakpkl' → use .peakpkl extension
                    # After splitext(), stem will be just "{mzML_stem}" matching PSM base_name
                    psms_out = str(library_dir / f"{original_stem}.psmpkl")
                    peaks_out = str(library_dir / f"{original_stem}.peakpkl")

                    convert_cmd = [
                        "easypqp", "convert",
                        "--pepxml", filter_idxml,
                        "--spectra", matching_mzml,
                        "--psms", psms_out,
                        "--peaks", peaks_out
                    ]
                    if self.executor.run_command(convert_cmd):
                        psms_files.append(psms_out)
                        peaks_files.append(peaks_out)

                if psms_files:
                    # easypqp library outputs TSV format (despite common .pqp extension)
                    library_tsv = str(library_dir / "spectral_library.tsv")
                    library_cmd = ["easypqp", "library", "--out", library_tsv]

                    if not self.params.get("library-use-fdr", False):
                        # --nofdr only skips FDR recalculation, NOT threshold filtering
                        # Set all thresholds to 1.0 to bypass filtering for pre-filtered input
                        library_cmd.extend([
                            "--nofdr",
                            "--psm_fdr_threshold", "1.0",
                            "--peptide_fdr_threshold", "1.0",
                            "--protein_fdr_threshold", "1.0"
                        ])
                    else:
                        # Apply user-specified FDR filtering
                        library_cmd.extend([
                            "--psm_fdr_threshold",
                            str(self.params.get("library-psm-fdr", 0.01)),
                            "--peptide_fdr_threshold",
                            str(self.params.get("library-peptide-fdr", 0.01)),
                            "--protein_fdr_threshold",
                            str(self.params.get("library-protein-fdr", 0.01))
                        ])

                    for psms, peaks in zip(psms_files, peaks_files):
                        library_cmd.extend([psms, peaks])

                    if self.executor.run_command(library_cmd):
                        self.logger.log("✅ Spectral library created")
                        st.success("Spectral library created")
                    else:
                        self.logger.log("Warning: Failed to build spectral library")
                else:
                    self.logger.log("Warning: No PSMs converted for library generation")

            st.success(f"✓ {stem} identification completed")

            # # ================================
            # # 4️⃣ ProteomicsLFQ (cross-sample)
            # # ================================
            self.logger.log("📈 Running cross-sample quantification...")
            st.info("Running ProteomicsLFQ (cross-sample quantification)")

            quant_mztab = str(quant_dir / "openms_quant.mzTab")
            quant_cxml = str(quant_dir / "openms.consensusXML")
            quant_msstats = str(quant_dir / "openms_msstats.csv")

            # A fractionated experiment needs the design, so ProteomicsLFQ combines
            # the fractions of each sample; otherwise every file is its own sample.
            design = load_design(self.workflow_dir, self.params["mzML-files"])
            # Converted vendor runs are named by their mzML in the design for ProteomicsLFQ
            names = {Path(a).name: Path(b).name for a, b in zip(self.params["mzML-files"], in_mzML)}
            design["file"] = design["file"].map(lambda f: names.get(f, f))
            lfq_design = {}
            if is_fractionated(design):
                design_file = save_design(quant_dir, design, in_mzML)
                lfq_design = {"design": str(design_file)}
                self.logger.log(f"Experimental design with fractions: {design_file}")

            if match_between_runs(self.params):
                self.logger.log("Match between runs: on (unidentified features are linked to identified ones)")

            with st.spinner("ProteomicsLFQ"):
                    combined_in = " ".join(in_mzML)
                    combined_ids = " ".join(filter_results)
                    self.logger.log(f"COMBINED_IN {combined_in}", 1)
                    self.logger.log(f"COMBINED_IN_TYPE {type(combined_in).__name__}", 1)
                    self.logger.log(f"FILTER_RESULTS = {filter_results}", 1)
                    self.logger.log(f"FILTER_RESULTS_LEN = {len(filter_results)}", 1)

                    # ✅ Streamlit output (debug view)
                    st.markdown("### 🔍 ProteomicsLFQ Input Debug")
                    st.write("**combined_in:**", combined_in)
                    st.write("**combined_in type:**", type(combined_in).__name__)

                    st.write("**combined_ids:**", combined_ids)
                    st.write("**combined_ids type:**", type(combined_ids).__name__)

                    if not self.executor.run_topp(
                            "ProteomicsLFQ",
                            {
                                "in": [in_mzML],
                                "ids": [filter_results],
                                "out": [quant_mztab],
                                "out_cxml": [quant_cxml],
                                "out_msstats": [quant_msstats],
                            },
                            {
                                "fasta": str(database_fasta),
                                **lfq_design,
                                **lfq_transfer_params(self.params),
                                "threads": 12,
                                # Disable FAIMS/IM handling to avoid segfault in OpenMS 3.5.0
                                "PeptideQuantification:extract:IM_window": "0.0",
                                "PeptideQuantification:faims:merge_features": "false",
                            },
                        ):
                        self.logger.log("Workflow stopped due to error")
                        return False
            self.logger.log("✅ Quantification complete")

            # if not Path(quant_mztab).exists():
            #     st.error("ProteomicsLFQ failed: mzTab not created")
            #     st.stop()

            self.write_qpx()
            self.write_identified_spectra(in_mzML, filter_results)

            # ================================
            # 5️⃣ Final report
            # # ================================
            st.success("🎉 TOPP workflow completed successfully")
            st.write("📁 Results directory:")   
            st.code(str(results_dir)) 
            
            st.write("📄 Generated files:")
            st.write(f"- mzTab: {quant_mztab}")
            st.write(f"- consensusXML: {quant_cxml}")
            st.write(f"- MSstats CSV: {quant_msstats}")

            return True
        else:
            self.logger.log("⚙️ Running TMT workflow")

            results_dir = Path(self.results_dir)
            iso_dir = results_dir / "isobaric_consensusXML"
            search_dir = results_dir / "search_results"
            perc_dir = results_dir / "percolator_results"
            psm_filter_dir = results_dir / "psm_filter"
            map_dir = results_dir / "idmapper"
            merge_dir = results_dir / "merged"
            protein_dir = results_dir / "protein"
            msstats_dir = results_dir / "msstats"
            quant_dir = results_dir / "quant_results"

            iso_consensus = []
            search_results = []
            percolator_results = []
            psm_filtered = []
            mapped_ids = []

            for d in [
                iso_dir, search_dir, perc_dir, psm_filter_dir,
                map_dir, merge_dir, protein_dir, msstats_dir, quant_dir
            ]:
                d.mkdir(parents=True, exist_ok=True)

            for mz in in_mzML:
                stem = Path(mz).stem
                iso_consensus.append(str(iso_dir / f"{stem}_iso.consensusXML"))
                search_results.append(str(search_dir / f"{stem}_prose.idXML"))
                percolator_results.append(str(perc_dir / f"{stem}_prose_perc.idXML"))
                psm_filtered.append(str(psm_filter_dir / f"{stem}_prose_perc_filter.idXML"))
                mapped_ids.append(str(map_dir / f"{stem}_prose_perc_filter_map.consensusXML"))

            merged_id = str(merge_dir / "ID_mapper_merge.consensusXML")
            protein_id = str(protein_dir / "ID_mapper_merge_epi.consensusXML")
            protein_filter = str(protein_dir / "ID_mapper_merge_epi_filter.consensusXML")
            protein_resolved = str(protein_dir / "ID_mapper_merge_epi_filter_resconf.consensusXML")
            consensus_out = str(quant_dir / "openms_design_protein_openms.csv") 

            # --- IsobaricAnalyzer ---
            self.logger.log("🏷️ Running isobaric labeling analysis...")
            with st.spinner("IsobaricAnalyzer"):
                if not self.executor.run_topp(
                    "IsobaricAnalyzer",
                    {
                        "in": in_mzML,
                        "out": iso_consensus,
                    },
                    tool_instance_name="IsobaricAnalyzer-TMT",
                ):
                    self.logger.log("Workflow stopped due to error")
                    return False
            self.logger.log("✅ IsobaricAnalyzer complete")

            # --- ProSE database search ---
            self.logger.log("🔎 Running peptide search...")
            with st.spinner("ProSE database search"):
                if not self.run_search("ProSE-TMT", in_mzML, search_results, database_fasta):
                    self.logger.log("Workflow stopped due to error")
                    return False
            self.logger.log("✅ ProSE search complete")

            # Fragment tolerance of the search, for the spectrum views
            frag_tol, frag_tol_is_ppm = self.search_tolerance("ProSE-TMT")

            if not self._has_identifications(search_results, "peptide search"):
                return False

            # Build visualization cache for the search results
            results_dir_path = Path(self.results_dir)
            cache_dir = results_dir_path / "insight_cache"
            cache_dir.mkdir(parents=True, exist_ok=True)

            # Spectra cache (once, shared by all stages): one parquet per input
            # mzML, read only for the idXML of that file, so memory follows one
            # run instead of the whole study. Every input mzML gets a fixed file
            # index, so each PSM links to its own run's spectra.
            filename_to_index = {Path(m).name: i for i, m in enumerate(in_mzML)}
            spectra_dir = Path(self.workflow_dir, "spectra-cache")
            shutil.rmtree(spectra_dir, ignore_errors=True)
            spectra_df = write_spectra_cache(spectra_dir, filename_to_index, in_mzML)

            for idxml_file in search_results:
                idxml_path = Path(idxml_file)

                # Table, Heatmap, SequenceView and LinePlot caches
                if not build_psm_views(
                    idxml_path, cache_dir, spectra_df, filename_to_index, frag_tol, frag_tol_is_ppm,
                    export_dir=None,
                ):
                    self.logger.log(f"⚠️ No identifications in {idxml_path.name}, skipping its result views")

            self.logger.log("✅ Peptide search complete")

            # --- Optional AlphaPeptDeep features ---
            percolator_in = self.add_prediction_features(search_results, in_mzML, "ProSE-TMT")
            if percolator_in is None:
                return False

            # --- PercolatorAdapter ---
            self.logger.log("📊 Running rescoring...")
            with st.spinner(f"PercolatorAdapter"):
                if not self.executor.run_topp(
                    "PercolatorAdapter",
                    {
                        "in": percolator_in,
                        "out": percolator_results,
                    },
                    {"decoy_pattern": decoy_string},  # Always propagated from upstream
                    tool_instance_name="PercolatorAdapter-TMT",
                ):
                    self.logger.log("Workflow stopped due to error")
                    return False
            if not self._has_identifications(percolator_results, "rescoring"):
                return False

            # Build visualization cache for Percolator results
            for idxml_file in percolator_results:
                idxml_path = Path(idxml_file)

                # Table, Heatmap, SequenceView and LinePlot caches
                if not build_psm_views(
                    idxml_path, cache_dir, spectra_df, filename_to_index, frag_tol, frag_tol_is_ppm,
                    export_dir=None,
                ):
                    self.logger.log(f"⚠️ No identifications in {idxml_path.name}, skipping its result views")

            self.logger.log("✅ PercolatorAdapter complete")
                
            # --- IDFilter ---
            self.logger.log("🔧 Filtering identifications...")
            with st.spinner(f"IDFilter"):
                if not self.filter_psms(percolator_results, psm_filtered, "IDFilter-strict"):
                    self.logger.log("Workflow stopped due to error")
                    return False
            self.logger.log("✅ IDFilter-strict complete")

            if not self._has_identifications(psm_filtered, "FDR filtering"):
                return False

            # Build visualization cache for Filter results
            for idxml_file in psm_filtered:
                idxml_path = Path(idxml_file)

                # Table, Heatmap, SequenceView and LinePlot caches and the TSV exports
                if not build_psm_views(
                    idxml_path, cache_dir, spectra_df, filename_to_index, frag_tol, frag_tol_is_ppm,
                    export_dir=results_dir_path / "exports",
                ):
                    self.logger.log(f"⚠️ No identifications in {idxml_path.name}, skipping its result views")
            shutil.rmtree(spectra_dir, ignore_errors=True)

            # --- IDMapper ---
            self.logger.log("🗺️ Mapping IDs to isobaric consensus features...")
            for iso, psm, mapped in zip(iso_consensus, psm_filtered, mapped_ids):
                iso_stem = Path(iso).stem
                with st.spinner(f"IDMapper ({iso_stem})"):
                    if not self.executor.run_topp(
                        "IDMapper",
                        {
                            "in": [iso],
                            "id": [psm],
                            "out": [mapped],
                        },
                        tool_instance_name="IDMapper-TMT",
                    ):
                        self.logger.log("Workflow stopped due to error")
                        return False
            self.logger.log("✅ IDMapper complete")

            # --- FileMerger ---
            self.logger.log("🔗 Merging mapped consensus files...")
            with st.spinner("FileMerger"):
                if not self.executor.run_topp(
                    "FileMerger",
                    {
                        "in": mapped_ids,
                        "out": [merged_id],
                    },
                    tool_instance_name="FileMerger-TMT",
                ):
                    self.logger.log("Workflow stopped due to error")
                    return False
            self.logger.log("✅ FileMerger complete")

            # --- ProteinInference ---
            self.logger.log("🧩 Running protein inference...")
            with st.spinner("ProteinInference"):
                if not self.executor.run_topp(
                    "ProteinInference",
                    {
                        "in": [merged_id],
                        "out": [protein_id],
                    },
                    {"picked_decoy_string": decoy_string},  # Always propagated from upstream
                    tool_instance_name="ProteinInference-TMT",
                ):
                    self.logger.log("Workflow stopped due to error")
                    return False
            self.logger.log("✅ ProteinInference complete")

            # --- IDFilter-lenient (Protein) ---    
            self.logger.log("🔬 Filtering proteins...")
            with st.spinner("IDFilter (Protein)"):
                if not self.executor.run_topp(
                    "IDFilter",
                    {
                        "in": [protein_id],
                        "out": [protein_filter],
                    },
                    tool_instance_name="IDFilter-lenient"
                ):
                    self.logger.log("Workflow stopped due to error")
                    return False
            self.logger.log("✅ IDFilter-lenient (Protein) complete")

            # ================================
            # ✨ NEW: 8️⃣ IDConflictResolver (protein_filter → protein_resolved)
            # ================================
            self.logger.log("⚖️ Resolving ID conflicts...")
            with st.spinner("IDConflictResolver"):
                if not self.executor.run_topp(
                    "IDConflictResolver",
                    {
                        "in": [protein_filter],
                        "out": [protein_resolved],
                    },
                    tool_instance_name="IDConflictResolver-TMT",
                ):
                    self.logger.log("Workflow stopped due to error")
                    return False
            self.logger.log("✅ IDConflictResolver complete")

            # ================================
            # ✨ NEW: 🔟 ProteinQuantifier (protein_resolved → consensus_out)
            # ================================
            self.logger.log("📐 Running protein quantification...")
            with st.spinner("ProteinQuantifier"):
                if not self.executor.run_topp(
                    "ProteinQuantifier",
                    {
                        "in": [protein_resolved],
                        "out": [consensus_out],
                    },
                    tool_instance_name="ProteinQuantifier-TMT",
                ):
                    self.logger.log("Workflow stopped due to error")
                    return False
            self.logger.log("✅ ProteinQuantifier complete")
            self.logger.log("📄 Generating protein table...")
            self.write_qpx()
            self.write_identified_spectra(in_mzML, psm_filtered)

            self.logger.log("🎉 WORKFLOW FINISHED")

    @st.fragment
    def results(self) -> None:

        st.title("📊 Results")

        comet_tab, perc_tab, filter_tab, lfq_tab = st.tabs([
            "🔍 Identification",
            "🔍 Rescoring",
            "🔍 Filtering",
            "🔍 Quantification"
        ])

        # ================================
        # 🔍 CometAdapter
        # ================================
        with comet_tab:

            comet_dir = search_results_dir(Path(self.results_dir))
            comet_files = sorted(comet_dir.glob("*.idXML"))

            if not comet_files:
                st.warning("⚠ No Identification output files found.")
                return

            selected_file = st.selectbox("📁 Select Identification result file", comet_files)

            def idxml_to_df(idxml_file):
                proteins, peptides = load_idxml(idxml_file)

                records = []
                for pep in peptides:
                    rt = pep.getRT()
                    mz = pep.getMZ()
                    for h in pep.getHits():
                        protein_refs = [ev.getProteinAccession() for ev in h.getPeptideEvidences()]
                        records.append({
                            "RT": rt,
                            "m/z": mz,
                            "Sequence": h.getSequence().toString(),
                            "Charge": h.getCharge(),
                            "Score": h.getScore(),
                            "Proteins": ",".join(protein_refs) if protein_refs else None,
                        })

                df = pd.DataFrame(records)
                if not df.empty:
                    df["Charge"] = df["Charge"].astype(str)
                    df["Charge_num"] = df["Charge"].astype(int)
                return df

            df = idxml_to_df(selected_file)

            if df.empty:
                st.info("No peptide hits found.")
                return
                
            st.dataframe(df, use_container_width=True)

            df_plot = df.reset_index()

            fig = px.scatter(
                df_plot,
                x="RT",
                y="m/z",
                color="Score",
                custom_data=["index", "Sequence", "Proteins"],
                color_continuous_scale=["#a6cee3", "#1f78b4", "#08519c", "#08306b"],
            )
            fig.update_traces(
                marker=dict(size=6, opacity=0.8),
                hovertemplate='<b>Index: %{customdata[0]}</b><br>'
                            + 'RT: %{x:.2f}<br>'
                            + 'm/z: %{y:.4f}<br>'
                            + 'Score: %{marker.color:.3f}<br>'
                            + 'Sequence: %{customdata[1]}<br>'
                            + 'Proteins: %{customdata[2]}<br>'
                            + '<extra></extra>'
            )
            fig.update_layout(
                coloraxis_colorbar=dict(title="Score"),
                hovermode="closest"
            )

            clicked = plotly_events(fig, click_event=True, hover_event=False, override_height=550, key="comet_plot")

            if clicked:
                row_index = clicked[0]["pointNumber"]
                st.subheader("📌 Selected Peptide Match")
                st.dataframe(df.iloc[[row_index]], use_container_width=True)

        # ================================
        # 🔍 PercolatorAdapter RESULTS
        # ================================
        with perc_tab:

            perc_dir = Path(self.results_dir, "percolator_results")
            perc_files = sorted(perc_dir.glob("*.idXML"))

            if not perc_files:
                st.warning("⚠ No Rescoring output files found.")
                return

            selected_perc = st.selectbox("📁 Select Rescoring result file", perc_files)

            df_p = idxml_to_df(selected_perc)

            if df_p.empty:
                st.info("No peptide hits found in Rescoring result.")
                return

            st.dataframe(df_p, use_container_width=True)

            df_plot_p = df_p.reset_index()

            fig2 = px.scatter(
                df_plot_p,
                x="RT",
                y="m/z",
                color="Score",
                custom_data=["index", "Sequence", "Proteins"],
                color_continuous_scale=["#a6cee3", "#1f78b4", "#08519c", "#08306b"],
            )
            fig2.update_traces(
                marker=dict(size=6, opacity=0.8),
                hovertemplate='<b>Index: %{customdata[0]}</b><br>'
                            + 'RT: %{x:.2f}<br>'
                            + 'm/z: %{y:.4f}<br>'
                            + 'Score: %{marker.color:.3f}<br>'
                            + 'Sequence: %{customdata[1]}<br>'
                            + 'Proteins: %{customdata[2]}<br>'
                            + '<extra></extra>'
            )
            fig2.update_layout(
                coloraxis_colorbar=dict(title="Score"),
                hovermode="closest"
            )

            clicked2 = plotly_events(fig2, click_event=True, hover_event=False, override_height=550, key="perc_plot")

            if clicked2:
                idx = clicked2[0]["pointNumber"]
                st.subheader("📌 Selected Rescoring Peptide Match")
                st.dataframe(df_p.iloc[[idx]], use_container_width=True)

        # ================================
        # 🔍 IDFilter RESULTS
        # ================================
        with filter_tab:

            filter_dir = Path(self.results_dir, "psm_filter")
            filter_files = sorted(filter_dir.glob("*.idXML"))

            if not filter_files:
                st.warning("⚠ No Filtering output files found.")
                return
            
            st.info("Here you can explore the PSM scatterplot along with the detailed PSM table.")

            selected_filter = st.selectbox("📁 Select Filtering result file", filter_files)

            df_f = idxml_to_df(selected_filter)

            if df_f.empty:
                st.info("No peptide hits found in Filtering result.")
                return

            st.dataframe(df_f, use_container_width=True)

            df_plot_f = df_f.reset_index()

            fig3 = px.scatter(
                df_plot_f,
                x="RT",
                y="m/z",
                color="Score",
                custom_data=["index", "Sequence", "Proteins"],
                color_continuous_scale=["#a6cee3", "#1f78b4", "#08519c", "#08306b"],
            )
            fig3.update_traces(
                marker=dict(size=6, opacity=0.8),
                hovertemplate='<b>Index: %{customdata[0]}</b><br>'
                            + 'RT: %{x:.2f}<br>'
                            + 'm/z: %{y:.4f}<br>'
                            + 'Score: %{marker.color:.3f}<br>'
                            + 'Sequence: %{customdata[1]}<br>'
                            + 'Proteins: %{customdata[2]}<br>'
                            + '<extra></extra>'
            )
            fig3.update_layout(
                coloraxis_colorbar=dict(title="Score"),
                hovermode="closest"
            )

            clicked3 = plotly_events(fig3, click_event=True, hover_event=False, override_height=550, key="filter_plot")

            if clicked3:
                idx3 = clicked3[0]["pointNumber"]
                st.subheader("📌 Selected Filtering Peptide Match")
                st.dataframe(df_f.iloc[[idx3]], use_container_width=True)

        # ================================
        # 📊 ProteomicsLFQ RESULTS
        # ================================
        with lfq_tab:

            results_dir = Path(self.results_dir)
            proteomicslfq_dir = results_dir / "quant_results"

            if not proteomicslfq_dir.exists():
                st.warning("❗ Quantification directory not found. Please run the analysis first.")
                return

            csv_files = sorted(proteomicslfq_dir.glob("*.csv"))

            if not csv_files:
                st.info("No CSV files found in the Quantification directory.")
                return

            csv_file = csv_files[0]

            # Protein / PSM table tab
            protein_tab, psm_tab, volcano_tab, pca_tab, heatmap_tab = st.tabs(["🧬 Protein Table", "📄 PSM-level Quantification Table", "🌋 Statistical Analysis Volcano Plot", "📊 PCA", "🔥 Heatmap",])

            try:
                df = pd.read_csv(csv_file)

                if df.empty:
                    st.info("No data found in this file.")
                    return

                # Protein-level Table
                with protein_tab:

                    st.markdown("### 🧬 Protein-Level Abundance Table")

                    st.info(
                        "💡INFO\n\n"
                        "This protein-level table is generated by grouping all PSMs that map to the "
                        "same protein and aggregating their intensities across samples.\n"
                        "Additionally, log2 fold change and p-values are calculated between sample groups."
                    )

                    # --------------------------------------------------
                    # 1️⃣ Load group information from self.params
                    # --------------------------------------------------
                    group_map = {
                        key[11:]: value  # Remove "mzML-group-" prefix
                        for key, value in self.params.items()
                        if key.startswith("mzML-group-") and value
                    }

                    if not group_map:
                        st.warning("No group information found. Please define sample groups first.")
                        return

                    # --------------------------------------------------
                    # 2️⃣ Sample → Group mapping
                    # --------------------------------------------------
                    df["Sample"] = df["Reference"].str.replace(".mzML", "", regex=False)
                    df["Group"] = df["Reference"].map(group_map)

                    df = df.dropna(subset=["Group"])

                    # --------------------------------------------------
                    # 3️⃣ Determine comparison groups
                    # --------------------------------------------------
                    groups = sorted(df["Group"].unique())

                    if len(groups) < 2:
                        # st.write("Group map:", group_map)
                        # st.write("Reference examples:", df["Reference"].unique()[:5])
                        # st.write("Sample examples:", df["Sample"].unique()[:5])
                        # st.write("Groups after mapping:", df["Group"].unique())
                        st.warning("At least two groups are required for statistical comparison.")
                        return

                    group1, group2 = groups[:2]

                    st.info(f"Statistical comparison: **{group2} vs {group1}**")

                    # --------------------------------------------------
                    # 4️⃣ Protein-level statistics
                    # --------------------------------------------------
                    stats_rows = []

                    for protein, protein_df in df.groupby("ProteinName"):
                        g1_vals = protein_df[protein_df["Group"] == group1]["Intensity"].values
                        g2_vals = protein_df[protein_df["Group"] == group2]["Intensity"].values

                        if len(g1_vals) < 2 or len(g2_vals) < 2:
                            pval = np.nan
                        else:
                            _, pval = ttest_ind(g1_vals, g2_vals, equal_var=False)

                        mean_g1 = np.mean(g1_vals) if len(g1_vals) > 0 else np.nan
                        mean_g2 = np.mean(g2_vals) if len(g2_vals) > 0 else np.nan

                        log2fc = np.log2(mean_g2 / mean_g1) if mean_g1 > 0 else np.nan

                        stats_rows.append(
                            {
                                "ProteinName": protein,
                                "log2FC": log2fc,
                                "p-value": pval,
                            }
                        )

                    stats_df = pd.DataFrame(stats_rows)

                    # --------------------------------------------------
                    # 5️⃣ Build protein abundance table (pivot)
                    # --------------------------------------------------
                    all_samples = sorted(df["Sample"].unique())
                    pivot_list = []

                    for protein, group_df in df.groupby("ProteinName"):
                        peptides = ";".join(group_df["PeptideSequence"].unique())
                        intensity_dict = group_df.groupby("Sample")["Intensity"].sum().to_dict()

                        intensity_dict_complete = {
                            sample: intensity_dict.get(sample, 0)
                            for sample in all_samples
                        }

                        row = {
                            "ProteinName": protein,
                            **intensity_dict_complete,
                            "PeptideSequence": peptides,
                        }

                        pivot_list.append(row)

                    pivot_df = pd.DataFrame(pivot_list)

                    # --------------------------------------------------
                    # 6️⃣ Merge statistics + reorder columns
                    # --------------------------------------------------
                    pivot_df = pivot_df.merge(stats_df, on="ProteinName", how="left")

                    pivot_df = pivot_df[
                        ["ProteinName", "log2FC", "p-value"] + all_samples + ["PeptideSequence"]
                    ]

                    # ================================
                    # Common matrix for PCA / Heatmap
                    # ================================
                    expr_df = pivot_df.set_index("ProteinName")[all_samples]
                    expr_df = expr_df.replace(0, np.nan)   # treat missing values (not detected)
                    expr_df = np.log2(expr_df + 1)         # log2 transform
                    expr_df = expr_df.dropna()             # remove proteins with missing values (row-wise)

                    st.session_state["pivot_df"] = pivot_df
                    st.session_state["expr_df"] = expr_df
                    st.session_state["group_map"] = group_map

                    # --------------------------------------------------
                    # 7️⃣ Display
                    # --------------------------------------------------
                    st.dataframe(
                        pivot_df.sort_values("p-value"),
                        use_container_width=True,
                    )

                # PSM-level Table
                with psm_tab:
                    st.markdown(f"### 📄 PSM-level Quantification Table")
                    st.info("💡INFO \n\n This table shows the PSM-level quantification data, including protein IDs,peptide sequences, charge states, and intensities across samples.Each row represents one peptide-spectrum match detected from the MS/MS analysis.")
                    st.dataframe(df, use_container_width=True)


            except Exception as e:
                st.error(f"Failed to load {csv_file.name}: {e}")