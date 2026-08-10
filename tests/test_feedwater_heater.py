import math

import pytest

pytest.importorskip("CoolProp")

from thermowave.components.feedwater_heater import FeedwaterHeater  # noqa: E402
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
    fwh = FeedwaterHeater(name="fwh")
    assert set(fwh.ports()) == {"hot_in", "hot_out", "cold_in", "cold_out"}
    assert fwh.report_category() == "phase_change"


def test_rejects_out_of_range_pr_and_subcool():
    with pytest.raises(ValueError, match="PR_hot"):
        FeedwaterHeater(name="fwh", PR_hot=0.0)
    with pytest.raises(ValueError, match="PR_cold"):
        FeedwaterHeater(name="fwh", PR_cold=1.5)
    with pytest.raises(ValueError, match="subcool"):
        FeedwaterHeater(name="fwh", subcool=-1.0)


def test_residuals_zero_and_energy_balanced_hand_calc():
    fwh = FeedwaterHeater(name="fwh", outlet_quality=0.0)
    # P_hot/mdot_hot/T_cold_in/mdot_cold chosen so the feedwater's outlet
    # temperature stays comfortably below the extraction steam's own
    # saturation temperature at P_hot (T_sat(1e5 Pa) ~ 372.8 K) -- a
    # physically feasible FWH spec, not an infeasible one (which would
    # legitimately show up as a negative reported pinch, exactly like
    # Condenser's own diagnostic pinch does for the same reason).
    P_hot, P_cold = 1.0e5, 1.0e6
    mdot_hot, mdot_cold = 1.0, 30.0

    h_hot_in = WATER.saturated_vapor_enthalpy(P_hot)  # extraction steam, sat vapor
    h_hot_out = WATER.saturated_liquid_enthalpy(P_hot)  # x=0 drain target
    Q = mdot_hot * (h_hot_in - h_hot_out)

    h_cold_in = WATER.enthalpy_pt(P_cold, 350.0)  # subcooled feedwater
    h_cold_out = h_cold_in + Q / mdot_cold

    state = _FakeState(
        fluid=WATER,
        node_values={
            "fwh.hot_in": (P_hot, h_hot_in), "fwh.hot_out": (P_hot, h_hot_out),
            "fwh.cold_in": (P_cold, h_cold_in), "fwh.cold_out": (P_cold, h_cold_out),
        },
        mdots={
            "fwh.hot_in": mdot_hot, "fwh.hot_out": mdot_hot,
            "fwh.cold_in": mdot_cold, "fwh.cold_out": mdot_cold,
        },
    )
    residuals = fwh.residuals(state)
    assert len(residuals) == 6
    for r in residuals:
        assert math.isclose(r, 0.0, abs_tol=1e-3)

    metrics = fwh.report_metrics(state)
    assert math.isclose(metrics["power [W]"], Q, rel_tol=1e-9)
    assert metrics["power [W]"] > 0.0  # heat given up by the extraction steam
    # Feedwater heat gain equals the hot-side duty -- ONE shared Q, not two
    # independently-fixed numbers with a disclosed gap between them (the
    # actual problem this component exists to fix).
    assert math.isclose(mdot_cold * (h_cold_out - h_cold_in), Q, rel_tol=1e-9)
    # TTD (pinch): saturation temp of the draining steam minus the
    # feedwater's own outlet temp -- must be positive for a physically
    # sane spec (feedwater can't leave hotter than the steam condensing on it).
    assert metrics["pinch [K]"] > 0.0


def test_subcool_mode_drain_below_saturation():
    fwh = FeedwaterHeater(name="fwh", subcool=5.0)
    P_hot = 3.0e5
    T_sat = WATER.saturation_temperature(P_hot)
    h_hot_out_expected = WATER.enthalpy_pt(P_hot, T_sat - 5.0)

    state = _FakeState(
        fluid=WATER,
        node_values={
            "fwh.hot_in": (P_hot, WATER.saturated_vapor_enthalpy(P_hot)),
            "fwh.hot_out": (P_hot, h_hot_out_expected),
            "fwh.cold_in": (1.0e6, WATER.enthalpy_pt(1.0e6, 400.0)),
            "fwh.cold_out": (1.0e6, WATER.enthalpy_pt(1.0e6, 400.0) + 1000.0),
        },
        mdots={"fwh.hot_in": 1.0, "fwh.hot_out": 1.0, "fwh.cold_in": 5.0, "fwh.cold_out": 5.0},
    )
    hot_energy_residual = fwh.residuals(state)[1]
    assert math.isclose(hot_energy_residual, 0.0, abs_tol=1e-6)


def test_end_to_end_two_stream_converges_with_shared_duty():
    hot_src = Source(name="hot_src", P=1.0e5, T=390.0, mdot=1.0)  # slightly superheated extraction
    cold_src = Source(name="cold_src", P=1.0e6, T=350.0, mdot=30.0)
    fwh = FeedwaterHeater(name="fwh", outlet_quality=0.0)
    hot_snk = Sink(name="hot_snk")
    cold_snk = Sink(name="cold_snk")

    network = Network(fluid=WATER)
    for c in (hot_src, cold_src, fwh, hot_snk, cold_snk):
        network.add_component(c)
    network.connect(hot_src, "out", fwh, "hot_in")
    network.connect(fwh, "hot_out", hot_snk, "in")
    network.connect(cold_src, "out", fwh, "cold_in")
    network.connect(fwh, "cold_out", cold_snk, "in")

    result = network.solve(tol=1e-6, max_iter=200)
    assert result.converged
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params, node_fluid=result.node_fluid,
    )
    metrics = fwh.report_metrics(state)

    P_hot_out, h_hot_out = result.node_P["fwh.hot_out"], result.node_h["fwh.hot_out"]
    assert math.isclose(h_hot_out, WATER.saturated_liquid_enthalpy(P_hot_out), rel_tol=1e-6)

    h_cold_in = result.node_h["fwh.cold_in"]
    h_cold_out = result.node_h["fwh.cold_out"]
    Q_cold = 30.0 * (h_cold_out - h_cold_in)
    assert math.isclose(Q_cold, metrics["power [W]"], rel_tol=1e-6)
    assert metrics["pinch [K]"] > 0.0
