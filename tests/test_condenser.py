import math

import pytest

pytest.importorskip("CoolProp")

from thermowave.components.condenser import Condenser  # noqa: E402
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


def test_ports_and_category():
    cd = Condenser(name="cd")
    assert set(cd.ports()) == {"wf_in", "wf_out", "cool_in", "cool_out"}
    assert cd.report_category() == "phase_change"


def test_residuals_zero_and_energy_balanced_hand_calc():
    cd = Condenser(name="cd", outlet_quality=0.0)
    P_wf, P_cool = 1.0e5, 2.0e5
    mdot_wf, mdot_cool = 1.0, 40.0

    h_wf_in = WATER.saturated_vapor_enthalpy(P_wf)  # incoming sat vapor
    h_wf_out = WATER.saturated_liquid_enthalpy(P_wf)  # x=0 target
    Q = mdot_wf * (h_wf_in - h_wf_out)

    h_cool_in = WATER.enthalpy_pt(P_cool, 300.0)  # cold coolant
    h_cool_out = h_cool_in + Q / mdot_cool

    state = _FakeState(
        fluid=WATER,
        node_values={
            "cd.wf_in": (P_wf, h_wf_in), "cd.wf_out": (P_wf, h_wf_out),
            "cd.cool_in": (P_cool, h_cool_in), "cd.cool_out": (P_cool, h_cool_out),
        },
        mdots={
            "cd.wf_in": mdot_wf, "cd.wf_out": mdot_wf,
            "cd.cool_in": mdot_cool, "cd.cool_out": mdot_cool,
        },
    )
    residuals = cd.residuals(state)
    assert len(residuals) == 6
    for r in residuals:
        assert math.isclose(r, 0.0, abs_tol=1e-3)

    metrics = cd.report_metrics(state)
    assert math.isclose(metrics["power [W]"], Q, rel_tol=1e-9)
    assert metrics["power [W]"] > 0.0  # heat rejected by wf
    # Coolant heat gain equals wf duty.
    assert math.isclose(mdot_cool * (h_cool_out - h_cool_in), Q, rel_tol=1e-9)


def test_end_to_end_two_stream_converges():
    wf_src = Source(name="wf_src", P=1.0e5, T=400.0, mdot=1.0)  # slightly superheated
    cool_src = Source(name="cool_src", P=2.0e5, T=300.0, mdot=40.0)
    cd = Condenser(name="cd", outlet_quality=0.0)
    wf_snk = Sink(name="wf_snk")
    cool_snk = Sink(name="cool_snk")

    network = Network(fluid=WATER)
    for c in (wf_src, cool_src, cd, wf_snk, cool_snk):
        network.add_component(c)
    network.connect(wf_src, "out", cd, "wf_in")
    network.connect(cd, "wf_out", wf_snk, "in")
    network.connect(cool_src, "out", cd, "cool_in")
    network.connect(cd, "cool_out", cool_snk, "in")

    result = network.solve(tol=1e-6, max_iter=200)
    assert result.converged
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params,
    )
    metrics = cd.report_metrics(state)
    assert math.isclose(metrics["x_out [-]"], 0.0, abs_tol=1e-4)
    assert metrics["pinch [K]"] > 0.0
