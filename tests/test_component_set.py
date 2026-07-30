import pytest

from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid


def _source() -> Source:
    return Source(name="src", P=1.0, T=300.0, mdot=1.0)


def test_set_updates_attribute_and_returns_self():
    src = _source()
    result = src.set(mdot=2.0)
    assert src.mdot == 2.0
    assert result is src


def test_set_accepts_multiple_kwargs():
    src = _source()
    src.set(mdot=2.0, mdot_guess=3.0)
    assert src.mdot == 2.0
    assert src.mdot_guess == 3.0


def test_set_rejects_unknown_parameter():
    src = _source()
    with pytest.raises(AttributeError):
        src.set(not_a_real_param=1.0)


def test_set_rejects_name():
    src = _source()
    with pytest.raises(AttributeError):
        src.set(name="renamed")


def test_set_rejects_private_attribute():
    src = _source()
    with pytest.raises(AttributeError):
        src.set(_outlet_node="whatever")


def test_set_before_adding_to_network_does_not_error():
    src = _source()
    src.set(mdot=2.0)  # no _network attribute yet; must be a no-op, not a crash


def test_set_after_add_component_invalidates_network_caches():
    net = Network(fluid=IdealGasFluid(name="air", R=287.05, cp=1005.0))
    src = _source()
    net.add_component(src)
    version_before = net._topology_version
    src.set(mdot=2.0)
    assert net._topology_version > version_before
