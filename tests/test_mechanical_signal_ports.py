"""Core plumbing for kind="mechanical"/kind="signal" connections (Network,
NetworkState, Solver) -- exercised directly with minimal dummy components
rather than any real turbomachinery, since this is infrastructure Part 1
(Turbine/Compressor/Shaft/Generator) builds on top of, not a test of any of
those components themselves.
"""

from thermowave.components.base_component import BaseComponent
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _FixedShaft(BaseComponent):
    """No flow ports at all -- just pins a mechanical port to a literal N."""

    def __init__(self, name: str, N: float):
        self.name = name
        self.N = N
        self._shaft = f"{name}.shaft"

    def ports(self) -> dict[str, str]:
        return {}

    def mechanical_ports(self) -> dict[str, str]:
        return {"shaft": self._shaft}

    def fixed_mechanical_values(self) -> dict[str, float]:
        return {self._shaft: self.N}

    def residuals(self, state: "NetworkState") -> list[float]:
        return []


class _FreeShaft(BaseComponent):
    """No flow ports -- leaves its mechanical port free (a Newton unknown
    unless connected to something that fixes it)."""

    def __init__(self, name: str, guess: float = 500.0):
        self.name = name
        self.guess = guess
        self._shaft = f"{name}.shaft"

    def ports(self) -> dict[str, str]:
        return {}

    def mechanical_ports(self) -> dict[str, str]:
        return {"shaft": self._shaft}

    def free_mechanical_ports(self) -> dict[str, float]:
        return {self._shaft: self.guess}

    def closes_mechanical_nodes(self) -> list[str]:
        # A real closing residual would tie this node to some target; empty
        # here on purpose so the "unclosed free mechanical node" wiring
        # check has something to catch in the dedicated test below.
        return []

    def residuals(self, state: "NetworkState") -> list[float]:
        return []


class _PowerProvider(BaseComponent):
    """No flow ports -- publishes a fixed value on a signal port."""

    def __init__(self, name: str, power: float):
        self.name = name
        self.power = power
        self._power_port = f"{name}.power"

    def ports(self) -> dict[str, str]:
        return {}

    def signal_ports(self) -> dict[str, str]:
        return {"power": self._power_port}

    def provided_signal_values(self, state: "NetworkState") -> dict[str, float]:
        return {self._power_port: self.power}

    def residuals(self, state: "NetworkState") -> list[float]:
        return []


class _PowerReader(BaseComponent):
    """No flow ports -- reads whatever its signal port resolves to and
    stashes it for the test to inspect (residuals() has no other use for
    it here)."""

    def __init__(self, name: str):
        self.name = name
        self._power_port = f"{name}.power"
        self.last_read: float | None = None

    def ports(self) -> dict[str, str]:
        return {}

    def signal_ports(self) -> dict[str, str]:
        return {"power": self._power_port}

    def residuals(self, state: "NetworkState") -> list[float]:
        self.last_read = state.signal(self._power_port)
        return []


def _flow_only_network() -> Network:
    """A trivial Source->Sink flow network, just to satisfy
    Network.validate_topology()'s "needs a boundary fixing (P, h)" check --
    the mechanical/signal components under test contribute no flow ports at
    all, so this is otherwise unrelated to what's being verified.
    """
    network = Network(fluid=AIR)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(snk)
    network.connect(src, "out", snk, "in")
    return network


def test_connecting_free_mechanical_port_to_fixed_one_needs_no_unknown():
    network = _flow_only_network()
    fixed = _FixedShaft("fx", N=3000.0)
    free = _FreeShaft("fr", guess=500.0)
    network.add_component(fixed)
    network.add_component(free)
    network.connect(fixed, "shaft", free, "shaft", kind="mechanical")

    result = network.solve(progress=False)

    assert result.converged
    state = result.state()
    assert state.N(fixed._shaft) == 3000.0
    assert state.N(free._shaft) == 3000.0


def test_unconnected_free_mechanical_port_is_a_newton_unknown_needing_closure():
    network = _flow_only_network()
    free = _FreeShaft("fr", guess=500.0)
    network.add_component(free)

    problems = network.check_wiring()
    assert any("mechanical node" in p and "nothing closes it" in p for p in problems)

    try:
        network.solve(progress=False)
        raised = False
    except Exception:
        raised = True
    assert raised  # under-determined: one free N, zero closing residuals


def test_signal_port_carries_a_computed_value_between_components():
    network = _flow_only_network()
    provider = _PowerProvider("gen", power=42000.0)
    reader = _PowerReader("load")
    network.add_component(provider)
    network.add_component(reader)
    network.connect(provider, "power", reader, "power", kind="signal")

    result = network.solve(progress=False)

    assert result.converged
    assert reader.last_read == 42000.0


def test_network_connect_mechanical_and_signal_share_the_flow_union_find_but_not_the_flow_graph():
    network = Network(fluid=AIR)
    fixed = _FixedShaft("fx", N=1000.0)
    free = _FreeShaft("fr")
    network.add_component(fixed)
    network.add_component(free)
    network.connect(fixed, "shaft", free, "shaft", kind="mechanical")

    # Mechanical ports never appear in the flow-node bookkeeping -- a
    # network with only mechanical/signal components and no Source has no
    # boundary (P, h) fixed at all.
    assert network._fixed_node_values() == {}
    assert network._all_nodes() == []
    assert network._fixed_mechanical_values() == {
        network._canonical(fixed._shaft): 1000.0
    }
