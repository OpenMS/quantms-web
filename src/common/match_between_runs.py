"""Match between runs for label-free quantification.

ProteomicsLFQ quantifies only identified peptides when ``targeted_only`` is
``true``. With ``false`` it also detects unidentified features in every run
and links them to identified ones after alignment, so a peptide identified
in one run is quantified in the others too (match between runs).

Streamlit-free: used by ``execution()`` in the queue worker.
"""

MATCH_BETWEEN_RUNS_KEY = "match-between-runs"


def match_between_runs(params: dict) -> bool:
    """True when the workflow's parameters switch match between runs on."""
    value = params.get(MATCH_BETWEEN_RUNS_KEY, False)
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return bool(value)


def lfq_transfer_params(params: dict) -> dict:
    """ProteomicsLFQ parameters for the chosen match-between-runs setting.

    Always explicit, so a ``targeted_only`` value stored from an earlier
    version of the configuration page cannot decide it.
    """
    return {"targeted_only": "false" if match_between_runs(params) else "true"}
