import math

import pytest

pytest.importorskip("CoolProp")

from thermowave.components.simple_evaporator import SimpleEvaporator  # noqa: E402
from thermowave.components.sink import Sink  # noqa: E402
from thermowave.components.source import Source  # noqa: E402
from thermowave.core.network import Network, NetworkState  # noqa: E402
from thermowave.fluids.ideal_gas import IdealGasFluid  # noqa: E402
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


def test_rejects_pr_out_of_range():
    with pytest.raises(ValueError, match="PR"):
        SimpleEvaporator(name="ev", PR=1.5)


def test_rejects_outlet_quality_out_of_range():
    with pytest.raises(ValueError, match="outlet_quality"):
        SimpleEvaporator(name="ev", outlet_quality=1.5)


def test_rejects_negative_superheat():
    with pytest.raises(ValueError, match="superheat"):
        SimpleEvaporator(name="ev", superheat=-5.0)


def test_ports_and_category():
    ev = SimpleEvaporator(name="ev")
    assert ev.ports() == {"in": "ev.in", "out": "ev.out"}
    assert ev.report_category() == "phase_change"


def test_requires_two_phase_fluid():
    ev = SimpleEvaporator(name="ev")
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    state = _FakeState(
        fluid=air,
        node_values={"ev.in": (1.0e5, 3.0e5), "ev.out": (1.0e5, 2.6e6)},
        mdots={"ev.in": 1.0, "ev.out": 1.0},
    )
    with pytest.raises(ValueError, match="two-phase-capable"):
        ev.residuals(state)


def test_residuals_zero_at_saturated_vapor_outlet_hand_calc():
    ev = SimpleEvaporator(name="ev", outlet_quality=1.0)
    P = 1.0e5
    h_in = WATER.enthalpy_pt(P, 350.0)
    h_out = WATER.saturated_vapor_enthalpy(P)  # x=1 target
    state = _FakeState(
        fluid=WATER,
        node_values={"ev.in": (P, h_in), "ev.out": (P, h_out)},
        mdots={"ev.in": 2.0, "ev.out": 2.0},
    )
    momentum, energy, mass = ev.residuals(state)
    assert math.isclose(momentum, 0.0, abs_tol=1e-6)
    assert math.isclose(energy, 0.0, abs_tol=1e-3)
    assert math.isclose(mass, 0.0, abs_tol=1e-12)


def test_superheat_target_matches_temperature_hand_calc():
    ev = SimpleEvaporator(name="ev", superheat=50.0)
    P = 1.0e5
    h_in = WATER.enthalpy_pt(P, 350.0)
    T_target = WATER.saturation_temperature(P) + 50.0
    h_out = WATER.enthalpy_pt(P, T_target)
    state = _FakeState(
        fluid=WATER,
        node_values={"ev.in": (P, h_in), "ev.out": (P, h_out)},
        mdots={"ev.in": 1.0, "ev.out": 1.0},
    )
    _momentum, energy, _mass = ev.residuals(state)
    assert math.isclose(energy, 0.0, abs_tol=1e-3)


def test_duty_mode_sets_outlet_from_heat_added():
    duty = 1.0e6  # W
    mdot = 2.0
    ev = SimpleEvaporator(name="ev", duty=duty)
    P = 1.0e5
    h_in = WATER.enthalpy_pt(P, 350.0)
    h_out = h_in + duty / mdot
    state = _FakeState(
        fluid=WATER,
        node_values={"ev.in": (P, h_in), "ev.out": (P, h_out)},
        mdots={"ev.in": mdot, "ev.out": mdot},
    )
    _momentum, energy, _mass = ev.residuals(state)
    assert math.isclose(energy, 0.0, abs_tol=1e-6)


def test_end_to_end_boils_water_to_saturated_vapor():
    src = Source(name="src", P=1.0e5, T=350.0, mdot=2.0)
    ev = SimpleEvaporator(name="ev", outlet_quality=1.0)
    snk = Sink(name="snk")

    network = Network(fluid=WATER)
    for c in (src, ev, snk):
        network.add_component(c)
    network.connect(src, "out", ev, "in")
    network.connect(ev, "out", snk, "in")

    result = network.solve(tol=1e-6, max_iter=100)
    assert result.converged
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params,
    )
    metrics = ev.report_metrics(state)
    assert math.isclose(metrics["x_out [-]"], 1.0, abs_tol=1e-4)
    assert metrics["power [W]"] > 0.0  # heat added
    # Duty should equal mdot * (h_g - h_in).
    expected_Q = 2.0 * (WATER.saturated_vapor_enthalpy(1.0e5) - WATER.enthalpy_pt(1.0e5, 350.0))
    assert math.isclose(metrics["power [W]"], expected_Q, rel_tol=1e-4)
