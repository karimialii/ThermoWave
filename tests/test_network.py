import pytest

from thermowave.components.pipe import Pipe
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.exceptions import NetworkTopologyError
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def _simple_network() -> Network:
    network = Network(fluid=AIR)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    pipe = Pipe(name="p1", L=5.0, D=0.2, f=0.02)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(pipe)
    network.add_component(snk)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")
    return network


def test_network_state_node_returns_p_and_h():
    state = NetworkState(
        fluid=AIR,
        node_P={"n1": 101325.0},
        node_h={"n1": 301500.0},
        node_mdot={"n1": 1.0},
    )
    assert state.node("n1") == (101325.0, 301500.0)
    assert state.mdot("n1") == 1.0


def test_network_connect_merges_ports_into_one_canonical_node():
    network = Network(fluid=AIR)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(snk)
    network.connect(src, "out", snk, "in")

    assert network._canonical("src.out") == network._canonical("snk.in")


def test_network_connect_rejects_unknown_port_name():
    network = Network(fluid=AIR)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(snk)
    with pytest.raises(NetworkTopologyError, match="no port"):
        network.connect(src, "wrong_port", snk, "in")


def test_network_connect_rejects_unsupported_kind():
    network = Network(fluid=AIR)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(snk)
    with pytest.raises(NotImplementedError, match="mechanical"):
        network.connect(src, "out", snk, "in", kind="mechanical")


def test_network_all_nodes_includes_ports_and_internal_nodes():
    network = Network(fluid=AIR)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    pipe = Pipe(name="p1", L=6.0, D=0.2, f=0.02, n_elem=3)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(pipe)
    network.add_component(snk)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")
    assert network._all_nodes() == ["src.out", "p1__mid1", "p1__mid2", "p1.out"]


def test_network_port_nodes_excludes_pipe_internal_nodes():
    network = Network(fluid=AIR)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    pipe = Pipe(name="p1", L=6.0, D=0.2, f=0.02, n_elem=3)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(pipe)
    network.add_component(snk)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")
    assert network._port_nodes() == ["src.out", "p1.out"]


def test_network_fixed_node_values_collected_from_source():
    network = _simple_network()
    fixed = network._fixed_node_values()
    assert set(fixed.keys()) == {"src.out"}
    P, h = fixed["src.out"]
    assert P == 101325.0
    assert h == 1005.0 * 300.0


def test_network_fixed_node_mdot_from_source():
    network = _simple_network()
    assert network._fixed_node_mdot() == {"src.out": 1.0}


def test_network_fixed_node_mdot_empty_when_no_source():
    network = Network(fluid=AIR)
    network.add_component(Sink(name="snk"))
    assert network._fixed_node_mdot() == {}


def test_network_validate_topology_allows_no_mdot_source():
    # A Source may legitimately leave mdot unfixed (mdot=None) so total mass
    # flow is solved for instead, closed by some other residual elsewhere
    # (e.g. a Sink pinning its inlet pressure) — validate_topology no longer
    # treats "nothing fixes mdot" as an error by itself.
    class _FixedValueOnlyComponent(Pipe):
        def fixed_node_values(self, fluid):
            h = fluid.enthalpy_pt(101325.0, 300.0)
            return {self._inlet_node: (101325.0, h)}

    network = Network(fluid=AIR)
    src = _FixedValueOnlyComponent(name="p1", L=5.0, D=0.2, f=0.02)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(snk)
    network.connect(src, "out", snk, "in")
    network.validate_topology()  # must not raise

    # The network is still genuinely under-constrained (nothing pins mdot at
    # all), so the solver's own square-system check must catch it instead.
    with pytest.raises(NetworkTopologyError, match="unknown"):
        network.solve()


def test_network_validate_topology_raises_when_no_boundary_state():
    network = Network(fluid=AIR)
    network.add_component(Pipe(name="p1", L=5.0, D=0.2, f=0.02))
    with pytest.raises(NetworkTopologyError, match="boundary"):
        network.validate_topology()


def test_network_validate_topology_passes_for_simple_network():
    network = _simple_network()
    network.validate_topology()  # must not raise


def test_source_mdot_none_closed_by_pressure_sink_solves_for_mdot():
    # A pipe pressure-drops with flow; pinning the exit to a lower ambient
    # pressure than the source, with mdot free, should solve for whatever
    # mdot produces exactly that drop across the pipe.
    network = Network(fluid=AIR)
    src = Source(name="src", P=200000.0, T=300.0, mdot=None, mdot_guess=0.5)
    pipe = Pipe(name="p1", L=5.0, D=0.2, f=0.02)
    snk = Sink(name="snk", P=150000.0)
    network.add_component(src)
    network.add_component(pipe)
    network.add_component(snk)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")

    result = network.solve(tol=1e-10, max_iter=200, damping=0.6)

    assert result.node_mdot["src.out"] > 0.0
    assert result.node_P["snk.in"] == pytest.approx(150000.0, rel=1e-6)
