"""Reporter channels of the isobaric label types IsobaricAnalyzer supports."""
import re

# Channel names in IsobaricAnalyzer's channel order, so channel n is the n-th entry
_TMT11 = ["126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C"]
_TMT16 = _TMT11 + ["132N", "132C", "133N", "133C", "134N"]
ISOBARIC_CHANNELS = {
    "itraq4plex": ["114", "115", "116", "117"],
    "itraq8plex": ["113", "114", "115", "116", "117", "118", "119", "121"],
    "tmt6plex": ["126", "127", "128", "129", "130", "131"],
    "tmt10plex": _TMT11[:9] + ["131"],
    "tmt11plex": _TMT11,
    "tmt16plex": _TMT16,
    "tmt18plex": _TMT16 + ["134C", "135N"],
}

_CHANNEL_COLUMN = re.compile(r"^(?:abundance_)?sample(\d+)\[[^\]]*\]$")


def label_type(params: dict) -> str:
    """The label type of the run, from the saved workflow parameters."""
    for source in (params.get("IsobaricAnalyzer-TMT", {}), params.get("_defaults", {}).get("IsobaricAnalyzer-TMT", {})):
        if source.get("type"):
            return str(source["type"]).lower()
    return "tmt11plex"


def channel_column_names(columns, label: str) -> dict:
    """Readable names for ProteinQuantifier's channel columns.

    ProteinQuantifier writes ``abundance_sample3[unknown]`` when the channels
    have no description; this becomes ``sample3[127C]``. The ``sampleN[`` part
    stays, because the group assignment finds channels by it.
    """
    tags = ISOBARIC_CHANNELS.get(label, [])
    renamed = {}
    for col in columns:
        m = _CHANNEL_COLUMN.match(str(col))
        if m:
            n = int(m.group(1))
            tag = tags[n - 1] if n <= len(tags) else str(col).rsplit("[", 1)[1].rstrip("]")
            renamed[col] = f"sample{n}[{tag}]"
    return renamed
