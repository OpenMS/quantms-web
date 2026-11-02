"""Add AlphaPeptDeep prediction features to search results before Percolator.

usage: python peptdeep_rescore.py <params.json>

params.json holds equal-length lists ``in`` (search idXMLs), ``mzml`` (their
spectra) and ``out`` (annotated idXMLs), plus ``fragment_tolerance``,
``fragment_ppm``, ``instrument`` ("auto" reads it from each mzML header) and
``summary`` (a JSON file the per-file results are written to).

The model is loaded once for all files, in this process only, so its memory is
returned when the step ends. See ``src/common/peptdeep_features.py``.
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

DEFAULTS = [
    {"key": "in", "value": [], "hide": True},
    {"key": "mzml", "value": [], "hide": True},
    {"key": "out", "value": [], "hide": True},
    {"key": "summary", "value": "", "hide": True},
    {"key": "fragment_tolerance", "value": 20.0, "hide": True},
    {"key": "fragment_ppm", "value": True, "hide": True},
    {"key": "instrument", "value": "auto", "hide": True},
]


def main(params_file: str) -> int:
    params = {entry["key"]: entry["value"] for entry in DEFAULTS}
    params.update(json.loads(Path(params_file).read_text(encoding="utf-8")))
    if not (len(params["in"]) == len(params["mzml"]) == len(params["out"])):
        print("in, mzml and out must have the same length", file=sys.stderr)
        return 2

    from src.common.ms2_prediction import load_model
    from src.common.peptdeep_features import add_features

    started = time.time()
    model = load_model()
    print(f"AlphaPeptDeep models loaded in {time.time() - started:.0f} s", flush=True)
    instrument = None if params["instrument"] in ("", "auto", None) else params["instrument"]
    summaries = []
    for idxml_in, mzml, idxml_out in zip(params["in"], params["mzml"], params["out"]):
        started = time.time()
        summary = add_features(
            idxml_in, mzml, idxml_out, model,
            tolerance=float(params["fragment_tolerance"]), ppm=bool(params["fragment_ppm"]),
            instrument=instrument, log=lambda line: print(line, flush=True),
        )
        summary["seconds"] = round(time.time() - started, 1)
        summaries.append(summary)
    if params["summary"]:
        Path(params["summary"]).write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
