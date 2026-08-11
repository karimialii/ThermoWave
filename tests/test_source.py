import pytest

from thermowave.components.source import Source


def test_source_rejects_nonpositive_pressure():
    try:
        Source(name="s1", P=0.0, T=300.0, mdot=1.0)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_source_rejects_nonpositive_temperature():
    try:
        Source(name="s1", P=1e5, T=0.0, mdot=1.0)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_source_rejects_negative_mdot():
    try:
        Source(name="s1", P=1e5, T=300.0, mdot=-1.0)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_source_accepts_valid_inputs():
    src = Source(name="s1", P=1e5, T=300.0, mdot=1.0)
    assert src.P_si == pytest.approx(1e5)
    assert src.mdot == 1.0


def test_source_accepts_none_mdot():
    src = Source(name="s1", P=1e5, T=300.0, mdot=None)
    assert src.mdot is None
