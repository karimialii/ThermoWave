import math

import pytest

pytest.importorskip("CoolProp")

from thermowave.components.evaporator import Evaporator  # noqa: E402
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
    ev = Evaporator(name="ev")
    assert set(ev.ports()) == {"wf_in", "wf_out", "src_in", "src_out"}
    assert ev.report_category() == "phase_change"


def test_residuals_zero_and_energy_balanced_hand_calc():
    ev = Evaporator(name="ev", outlet_quality=1.0)
    P_wf, P_src = 1.0e5, 5.0e5
    mdot_wf, mdot_src = 1.0, 30.0

    h_wf_in = WATER.enthalpy_pt(P_wf, 350.0)
    h_wf_out = WATER.saturated_vapor_enthalpy(P_wf)
    Q = mdot_wf * (h_wf_out - h_wf_in)

    h_src_in = WATER.enthalpy_pt(P_src, 420.0)
    h_src_out = h_src_in - Q / mdot_src

    state = _FakeState(
        fluid=WATER,
        node_values={
            "ev.wf_in": (P_wf, h_wf_in), "ev.wf_out": (P_wf, h_wf_out),
            "ev.src_in": (P_src, h_src_in), "ev.src_out": (P_src, h_src_out),
        },
        mdots={
            "ev.wf_in": mdot_wf, "ev.wf_out": mdot_wf,
            "ev.src_in": mdot_src, "ev.src_out": mdot_src,
        },
    )
    residuals = ev.residuals(state)
    assert len(residuals) == 6
    for r in residuals:
        assert math.isclose(r, 0.0, abs_tol=1e-3)

    metrics = ev.report_metrics(state)
    assert math.isclose(metrics["power [W]"], Q, rel_tol=1e-9)
    # Source-side heat loss equals wf duty (energy balance).
    assert math.isclose(mdot_src * (h_src_in - h_src_out), Q, rel_tol=1e-9)


def test_pinch_is_positive_for_a_feasible_spec():
    ev = Evaporator(name="ev", outlet_quality=1.0)
    P_wf, P_src = 1.0e5, 5.0e5
    mdot_wf, mdot_src = 1.0, 30.0
    h_wf_in = WATER.enthalpy_pt(P_wf, 350.0)
    h_wf_out = WATER.saturated_vapor_enthalpy(P_wf)
    Q = mdot_wf * (h_wf_out - h_wf_in)
    h_src_in = WATER.enthalpy_pt(P_src, 420.0)
    h_src_out = h_src_in - Q / mdot_src
    state = _FakeState(
        fluid=WATER,
        node_values={
            "ev.wf_in": (P_wf, h_wf_in), "ev.wf_out": (P_wf, h_wf_out),
            "ev.src_in": (P_src, h_src_in), "ev.src_out": (P_src, h_src_out),
        },
        mdots={
            "ev.wf_in": mdot_wf, "ev.wf_out": mdot_wf,
            "ev.src_in": mdot_src, "ev.src_out": mdot_src,
        },
    )
    assert ev.report_metrics(state)["pinch [K]"] > 0.0


def test_end_to_end_two_stream_converges_and_balances():
    wf_src = Source(name="wf_src", P=1.0e5, T=350.0, mdot=1.0)
    heat_src = Source(name="heat_src", P=5.0e5, T=420.0, mdot=30.0)
    ev = Evaporator(name="ev", outlet_quality=1.0)
    wf_snk = Sink(name="wf_snk")
    src_snk = Sink(name="src_snk")

    network = Network(fluid=WATER)
    for c in (wf_src, heat_src, ev, wf_snk, src_snk):
        network.add_component(c)
    network.connect(wf_src, "out", ev, "wf_in")
    network.connect(ev, "wf_out", wf_snk, "in")
    network.connect(heat_src, "out", ev, "src_in")
    network.connect(ev, "src_out", src_snk, "in")

    result = network.solve(tol=1e-6, max_iter=200)
    assert result.converged
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params,
    )
    metrics = ev.report_metrics(state)
    assert math.isclose(metrics["x_out [-]"], 1.0, abs_tol=1e-4)
    assert metrics["power [W]"] > 0.0

    # Energy balance across the exchanger.
    Q_src = 30.0 * (
        result.node_h["ev.src_in"] - result.node_h["ev.src_out"]
    )
    assert math.isclose(Q_src, metrics["power [W]"], rel_tol=1e-5)
