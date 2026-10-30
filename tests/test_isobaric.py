"""Tests for the reporter-channel naming in src/common/isobaric.py."""
from src.common.isobaric import ISOBARIC_CHANNELS, channel_column_names, label_type


def test_plex_sizes_match_their_names():
    for label, channels in ISOBARIC_CHANNELS.items():
        assert len(channels) == int(label.split("plex")[0].lstrip("tmitraq"))
        assert len(set(channels)) == len(channels)


def test_protein_quantifier_columns_get_reporter_names():
    columns = ["protein", "n_proteins", "abundance_sample1[unknown]", "abundance_sample3[unknown]"]
    assert channel_column_names(columns, "tmt11plex") == {
        "abundance_sample1[unknown]": "sample1[126]",
        "abundance_sample3[unknown]": "sample3[127C]",
    }


def test_described_channels_keep_their_description_beyond_the_plex():
    assert channel_column_names(["abundance_sample5[pool]"], "itraq4plex") == {
        "abundance_sample5[pool]": "sample5[pool]"
    }


def test_label_type_prefers_the_saved_choice():
    assert label_type({}) == "tmt11plex"
    assert label_type({"_defaults": {"IsobaricAnalyzer-TMT": {"type": "tmt16plex"}}}) == "tmt16plex"
    assert label_type({"IsobaricAnalyzer-TMT": {"type": "TMT6plex"}}) == "tmt6plex"
