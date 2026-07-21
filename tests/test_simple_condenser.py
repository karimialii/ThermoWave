import math

import pytest

pytest.importorskip("CoolProp")

from thermowave.components.simple_condenser import SimpleCondenser  # noqa: E402
from thermowave.components.sink import Sink  # noqa: E402
from thermowave.components.source import Source  # noqa: E402
from thermowave.core.network import Network, NetworkState  # noqa: E402
from thermowave.fluids.real_fluid import CoolPropFluid  # noqa: E402

WATER = CoolPropFluid(name="Water")


class _FakeState:
    def __init__(self, fluid, node_values, mdots):
        self.fluid = fluid
        self._node_values = node_values
        self._mdots = mdots

    def node(self, name):
        return self._node_values[name]

    def fluid_at(self, name):
        return self.fluid

    def mdot(self, name):
        return self._mdots[name]


def test_rejects_subcool_negative():
    with pytest.raises(ValueError, match="subcool"):
        SimpleCondenser(name="cd", subcool=-1.0)


def test_ports_and_category():
    cd = SimpleCondenser(name="cd")
    assert cd.ports() == {"in": "cd.in", "out": "cd.out"}
    assert cd.report_category() == "phase_change"


def test_residuals_zero_at_saturated_liquid_outlet_hand_calc():
    cd = SimpleCondenser(name="cd", outlet_quality=0.0)
    P = 1.0e5
    h_in = WATER.saturated_vapor_enthalpy(P)  # incoming saturated vapor
    h_out = WATER.saturated_liquid_enthalpy(P)  # x=0 target
    state = _FakeState(
        fluid=WATER,
        node_values={"cd.in": (P, h_in), "cd.out": (P, h_out)},
        mdots={"cd.in": 3.0, "cd.out": 3.0},
    )
    momentum, energy, mass = cd.residuals(state)
    assert math.isclose(momentum, 0.0, abs_tol=1e-6)
    assert math.isclose(energy, 0.0, abs_tol=1e-3)
    assert math.isclose(mass, 0.0, abs_tol=1e-12)


def test_power_is_negative_heat_rejected():
    cd = SimpleCondenser(name="cd", outlet_quality=0.0)
    P = 1.0e5
    h_in = WATER.saturated_vapor_enthalpy(P)
    h_out = WATER.saturated_liquid_enthalpy(P)
    state = _FakeState(
        fluid=WATER,
        node_values={"cd.in": (P, h_in), "cd.out": (P, h_out)},
        mdots={"cd.in": 3.0, "cd.out": 3.0},
    )
    metrics = cd.report_metrics(state)
    assert metrics["power [W]"] < 0.0
    assert math.isclose(metrics["x_out [-]"], 0.0, abs_tol=1e-4)


def test_subcool_target_is_below_saturation():
    cd = SimpleCondenser(name="cd", subcool=10.0)
    P = 1.0e5
    T_target = WATER.saturation_temperature(P) - 10.0
    h_out = WATER.enthalpy_pt(P, T_target)
    h_in = WATER.saturated_vapor_enthalpy(P)
    state = _FakeState(
        fluid=WATER,
        node_values={"cd.in": (P, h_in), "cd.out": (P, h_out)},
        mdots={"cd.in": 1.0, "cd.out": 1.0},
    )
    _momentum, energy, _mass = cd.residuals(state)
    assert math.isclose(energy, 0.0, abs_tol=1e-3)


def test_end_to_end_condenses_steam_to_saturated_liquid():
    # Feed slightly superheated steam in; condense to saturated liquid.
    src = Source(name="src", P=1.0e5, T=450.0, mdot=1.5)
    cd = SimpleCondenser(name="cd", outlet_quality=0.0)
    snk = Sink(name="snk")

    network = Network(fluid=WATER)
    for c in (src, cd, snk):
        network.add_component(c)
    network.connect(src, "out", cd, "in")
    network.connect(cd, "out", snk, "in")

    result = network.solve(tol=1e-6, max_iter=100)
    assert result.converged
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params,
    )
    metrics = cd.report_metrics(state)
    assert math.isclose(metrics["x_out [-]"], 0.0, abs_tol=1e-4)
    assert metrics["power [W]"] < 0.0
