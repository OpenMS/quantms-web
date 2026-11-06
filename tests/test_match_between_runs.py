from src.common.match_between_runs import MATCH_BETWEEN_RUNS_KEY, lfq_transfer_params, match_between_runs


def test_off_by_default():
    assert not match_between_runs({})
    assert lfq_transfer_params({}) == {"targeted_only": "true"}


def test_on():
    params = {MATCH_BETWEEN_RUNS_KEY: True}
    assert match_between_runs(params)
    assert lfq_transfer_params(params) == {"targeted_only": "false"}


def test_string_values():
    assert match_between_runs({MATCH_BETWEEN_RUNS_KEY: "true"})
    assert not match_between_runs({MATCH_BETWEEN_RUNS_KEY: "false"})


def test_configure_page_no_longer_pins_targeted_only():
    """The checkbox decides; a fixed custom default would hide that."""
    from pathlib import Path

    source = Path("src/WorkflowTest.py").read_text()
    assert '"targeted_only": "true"' not in source
    assert "**lfq_transfer_params(self.params)" in source
