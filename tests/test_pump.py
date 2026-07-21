import math

import pytest

pytest.importorskip("CoolProp")

from thermowave.components.pump import Pump  # noqa: E402
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


def test_rejects_both_or_neither_of_p_out_and_pr():
    with pytest.raises(ValueError, match="exactly one"):
        Pump(name="p", P_out=1.0e6, PR=10.0)
    with pytest.raises(ValueError, match="exactly one"):
        Pump(name="p")


def test_rejects_pr_not_greater_than_one():
    with pytest.raises(ValueError, match="PR"):
        Pump(name="p", PR=0.5)


def test_requires_entropy_fluid():
    p = Pump(name="p", P_out=1.0e6)
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    state = _FakeState(
        fluid=air,
        node_values={"p.in": (1.0e4, 2.0e5), "p.out": (1.0e6, 2.0e5)},
        mdots={"p.in": 1.0, "p.out": 1.0},
    )
    with pytest.raises(ValueError, match="entropy"):
        p.residuals(state)


def test_isentropic_rise_hand_calc_with_efficiency():
    eta = 0.75
    p = Pump(name="p", P_out=1.0e6, eta=eta)
    P_in = 1.0e4
    h_in = WATER.saturated_liquid_enthalpy(P_in)
    P_out = 1.0e6
    s_in = WATER.entropy_ph(P_in, h_in)
    h_out_isentropic = WATER.enthalpy_ps(P_out, s_in)
    h_out = h_in + (h_out_isentropic - h_in) / eta
    state = _FakeState(
        fluid=WATER,
        node_values={"p.in": (P_in, h_in), "p.out": (P_out, h_out)},
        mdots={"p.in": 1.0, "p.out": 1.0},
    )
    momentum, energy, mass = p.residuals(state)
    assert math.isclose(momentum, 0.0, abs_tol=1e-3)
    assert math.isclose(energy, 0.0, abs_tol=1e-6)
    assert math.isclose(mass, 0.0, abs_tol=1e-12)
    # Actual rise exceeds ideal for eta < 1.
    assert (h_out - h_in) > (h_out_isentropic - h_in)


def test_end_to_end_pump_raises_pressure_with_small_work():
    src = Source(name="src", P=1.0e4, T=300.0, mdot=2.0)  # subcooled liquid
    pump = Pump(name="pump", P_out=1.0e6, eta=0.75)
    snk = Sink(name="snk")

    network = Network(fluid=WATER)
    for c in (src, pump, snk):
        network.add_component(c)
    network.connect(src, "out", pump, "in")
    network.connect(pump, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=100)
    assert result.converged
    assert math.isclose(result.node_P["pump.out"], 1.0e6, rel_tol=1e-9)
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params,
    )
    metrics = pump.report_metrics(state)
    assert metrics["power [W]"] > 0.0  # work input
    # Pumping liquid is cheap: specific work << turbine-scale enthalpies.
    assert metrics["power [W]"] / 2.0 < 5.0e3  # < 5 kJ/kg
