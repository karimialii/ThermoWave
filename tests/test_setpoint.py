import math
from pathlib import Path

import pytest

from thermowave.components.compressor import Compressor
from thermowave.components.setpoint import Setpoint
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.exceptions import NetworkTopologyError
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

_MAP_PATH = str(Path(__file__).parent / "fixtures" / "simple_compressor_map.cop")

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def test_setpoint_raises_when_target_component_has_no_matching_free_parameter():
    comp = Compressor(name="c1", map_path=_MAP_PATH, N=3000.0, gamma=1.4)  # N fixed, not free
    with pytest.raises(ValueError, match="doesn't currently declare"):
        Setpoint(
            name="sp1", component=comp, free_param="N", target_metric="PR [-]", value=3.0
        )


def test_setpoint_ports_is_empty():
    comp = Compressor(name="c1", map_path=_MAP_PATH, gamma=1.4)  # N left free
    sp = Setpoint(name="sp1", component=comp, free_param="N", target_metric="PR [-]", value=3.0)
    assert sp.ports() == {}


def test_setpoint_residuals_is_metric_minus_value():
    comp = Compressor(name="c1", map_path=_MAP_PATH, gamma=1.4)
    sp = Setpoint(
        name="sp1", component=comp, free_param="N", target_metric="power [W]", value=500.0
    )

    class _FakeComponentState:
        pass

    class _StubComponent:
        name = "c1"

        def report_metrics(self, state):
            return {"power [W]": 700.0}

    sp.component = _StubComponent()
    assert sp.residuals(_FakeComponentState()) == [200.0]


def test_setpoint_drives_compressor_pr_target_end_to_end():
    gamma = 1005.0 / (1005.0 - 287.05)
    src = Source(name="src", P=101325.0, T=300.0, mdot=2.0 * 1.01325 / math.sqrt(300.0))
    comp = Compressor(name="c1", map_path=_MAP_PATH, gamma=gamma)  # N left free
    sp = Setpoint(
        name="sp1", component=comp, free_param="N", target_metric="PR [-]", value=3.0
    )
    snk = Sink(name="snk")

    network = Network(fluid=AIR)
    for component in (src, comp, sp, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=100, verbose=False)
    assert result.converged
    state = NetworkState(
        fluid=result.fluid,
        node_P=result.node_P,
        node_h=result.node_h,
        node_mdot=result.node_mdot,
        params=result.params,
    )
    metrics = comp.report_metrics(state)
    assert math.isclose(metrics["PR [-]"], 3.0, rel_tol=1e-4)


def test_network_raises_clear_error_when_underdetermined():
    gamma = 1005.0 / (1005.0 - 287.05)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    comp = Compressor(name="c1", map_path=_MAP_PATH, gamma=gamma)  # N left free, no Setpoint
    snk = Sink(name="snk")

    network = Network(fluid=AIR)
    for component in (src, comp, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", snk, "in")

    with pytest.raises(NetworkTopologyError, match="not solvable"):
        network.solve(tol=1e-8, max_iter=100, verbose=False)


def test_network_solve_verbose_prints_system_summary(capsys):
    gamma = 1005.0 / (1005.0 - 287.05)
    src = Source(name="src", P=101325.0, T=300.0, mdot=2.0 * 1.01325 / math.sqrt(300.0))
    comp = Compressor(name="c1", map_path=_MAP_PATH, gamma=gamma)
    sp = Setpoint(
        name="sp1", component=comp, free_param="N", target_metric="PR [-]", value=3.0
    )
    snk = Sink(name="snk")

    network = Network(fluid=AIR)
    for component in (src, comp, sp, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", snk, "in")

    network.solve(tol=1e-8, max_iter=100, verbose=True)
    out = capsys.readouterr().out
    assert "Network system:" in out
    assert "square, solvable" in out
