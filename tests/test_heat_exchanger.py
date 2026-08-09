import math

import pytest

from thermowave.components.heat_exchanger import (
    _CP_SAT_NUDGE_K,
    HeatExchanger,
    MultiPassHeatExchanger,
    _inlet_cp,
)
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _FakeState:
    def __init__(self, fluid, mdot: dict[str, float], node_values: dict[str, tuple[float, float]]):
        self.fluid = fluid
        self._mdot = mdot
        self._node_values = node_values

    def node(self, name: str) -> tuple[float, float]:
        return self._node_values[name]

    def fluid_at(self, name: str):
        return self.fluid

    def mdot(self, name: str) -> float:
        return self._mdot[name]


def _balanced_effectiveness_case(effectiveness=0.7, PR_hot=0.95, PR_cold=0.9):
    mdot_hot, mdot_cold = 1.0, 1.2
    P_hot, T_hot_in = 300000.0, 500.0
    P_cold, T_cold_in = 300000.0, 300.0
    h_hot_in = AIR.enthalpy_pt(P_hot, T_hot_in)
    h_cold_in = AIR.enthalpy_pt(P_cold, T_cold_in)

    hx = HeatExchanger(
        name="hx1", PR_hot=PR_hot, PR_cold=PR_cold, effectiveness=effectiveness
    )

    C_min = min(mdot_hot, mdot_cold) * 1005.0
    Q = effectiveness * C_min * (T_hot_in - T_cold_in)

    h_hot_out = h_hot_in - Q / mdot_hot
    h_cold_out = h_cold_in + Q / mdot_cold
    P_hot_out = PR_hot * P_hot
    P_cold_out = PR_cold * P_cold

    state = _FakeState(
        fluid=AIR,
        mdot={
            "hx1.hot_in": mdot_hot, "hx1.hot_out": mdot_hot,
            "hx1.cold_in": mdot_cold, "hx1.cold_out": mdot_cold,
        },
        node_values={
            "hx1.hot_in": (P_hot, h_hot_in), "hx1.hot_out": (P_hot_out, h_hot_out),
            "hx1.cold_in": (P_cold, h_cold_in), "hx1.cold_out": (P_cold_out, h_cold_out),
        },
    )
    return hx, state, Q


# --- _inlet_cp: saturation-boundary clamp -----------------------------------


def test_inlet_cp_passes_through_unclamped_for_non_two_phase_fluid():
    # AIR has no saturation dome (supports_two_phase() is False) -- nothing
    # to clamp away from, so _inlet_cp must be a plain pass-through to cp().
    P, T = 300000.0, 500.0
    assert _inlet_cp(AIR, P, T) == AIR.cp(P, T)


def test_inlet_cp_clamps_to_liquid_side_just_below_saturation():
    pytest.importorskip("CoolProp")
    from thermowave.fluids.real_fluid import CoolPropFluid

    water = CoolPropFluid(name="Water")
    P = 1.0e6
    T_sat = water.saturation_temperature(P)

    # Arbitrarily close to T_sat from below -- must land on the same clamped
    # liquid-side value regardless of how close, not flip based on whether
    # the raw cp() call happened to throw.
    for offset in (1e-2, 1e-4, 1e-8):
        T = T_sat - offset
        clamped = T_sat - _CP_SAT_NUDGE_K
        assert math.isclose(_inlet_cp(water, P, T), water.cp(P, clamped), rel_tol=1e-12)


def test_inlet_cp_clamps_to_vapor_side_just_above_saturation():
    pytest.importorskip("CoolProp")
    from thermowave.fluids.real_fluid import CoolPropFluid

    water = CoolPropFluid(name="Water")
    P = 1.0e6
    T_sat = water.saturation_temperature(P)

    for offset in (1e-2, 1e-4, 1e-8):
        T = T_sat + offset
        clamped = T_sat + _CP_SAT_NUDGE_K
        assert math.isclose(_inlet_cp(water, P, T), water.cp(P, clamped), rel_tol=1e-12)


def test_inlet_cp_unclamped_far_from_saturation_either_side():
    pytest.importorskip("CoolProp")
    from thermowave.fluids.real_fluid import CoolPropFluid

    water = CoolPropFluid(name="Water")
    P = 1.0e6
    T_sat = water.saturation_temperature(P)

    T_subcooled = T_sat - 50.0
    T_superheated = T_sat + 50.0
    assert math.isclose(_inlet_cp(water, P, T_subcooled), water.cp(P, T_subcooled), rel_tol=1e-9)
    assert math.isclose(
        _inlet_cp(water, P, T_superheated), water.cp(P, T_superheated), rel_tol=1e-9
    )
    # And the two sides genuinely differ (liquid vs. vapor cp) -- the clamp
    # doesn't paper over a real physical discontinuity, only the artificial
    # one right at the boundary.
    assert not math.isclose(
        _inlet_cp(water, P, T_subcooled), _inlet_cp(water, P, T_superheated), rel_tol=0.05
    )


def test_inlet_cp_uses_h_to_resolve_the_exactly_saturated_tie():
    pytest.importorskip("CoolProp")
    from thermowave.fluids.real_fluid import CoolPropFluid

    water = CoolPropFluid(name="Water")
    P = 1.0e6
    T_sat = water.saturation_temperature(P)
    h_f = water.saturated_liquid_enthalpy(P)
    h_g = water.saturated_vapor_enthalpy(P)

    # Both states sit at exactly T == T_sat -- indistinguishable in (P, T) --
    # so only h can say which side each is on. Saturated liquid must resolve
    # to the liquid-side cp, saturated vapor to the vapor-side cp.
    cp_liquid = water.cp(P, T_sat - _CP_SAT_NUDGE_K)
    cp_vapor = water.cp(P, T_sat + _CP_SAT_NUDGE_K)
    assert math.isclose(_inlet_cp(water, P, T_sat, h_f), cp_liquid, rel_tol=1e-12)
    assert math.isclose(_inlet_cp(water, P, T_sat, h_g), cp_vapor, rel_tol=1e-12)
    # ...and they are genuinely different sides, not both collapsed onto one.
    assert not math.isclose(cp_liquid, cp_vapor, rel_tol=0.05)


# --- construction / mode selection ----------------------------------------


def test_rejects_both_effectiveness_and_ua():
    with pytest.raises(ValueError, match="at most one"):
        HeatExchanger(name="hx1", PR_hot=1.0, PR_cold=1.0, effectiveness=0.7, UA=500.0)


def test_neither_effectiveness_nor_ua_means_free_ua():
    # UA mode with UA left as a free Newton unknown (see HeatExchanger's own
    # docstring) -- no longer an error, the same free-parameter pattern
    # SimpleHeatExchanger's own Q=None already uses.
    hx = HeatExchanger(name="hx1", PR_hot=1.0, PR_cold=1.0)
    assert hx._mode == "ua"
    assert hx.free_parameters() == {"UA": hx.UA_guess}


def test_rejects_arrangement_kwargs_in_effectiveness_mode():
    with pytest.raises(ValueError, match="only apply in UA mode"):
        HeatExchanger(
            name="hx1", PR_hot=1.0, PR_cold=1.0, effectiveness=0.7, n_passes=2
        )


def test_ports_returns_four_ports_in_either_mode():
    expected = {
        "hot_in": "hx1.hot_in", "hot_out": "hx1.hot_out",
        "cold_in": "hx1.cold_in", "cold_out": "hx1.cold_out",
    }
    assert HeatExchanger(name="hx1", PR_hot=0.95, PR_cold=0.9, effectiveness=0.7).ports() == expected
    assert HeatExchanger(name="hx1", PR_hot=0.95, PR_cold=0.9, UA=500.0).ports() == expected


# --- effectiveness mode ----------------------------------------------------


def test_effectiveness_mode_residuals_all_zero_at_exact_hand_calc_solution():
    hx, state, _Q = _balanced_effectiveness_case()
    residuals = hx.residuals(state)
    assert len(residuals) == 6
    for r in residuals:
        assert math.isclose(r, 0.0, abs_tol=1e-6)


def test_effectiveness_mode_duty_scales_linearly_with_effectiveness():
    hx_half, state_half, _ = _balanced_effectiveness_case(effectiveness=0.5)
    hx_full, state_full, _ = _balanced_effectiveness_case(effectiveness=1.0)
    assert math.isclose(
        hx_full._duty_fixed_effectiveness(state_full),
        2.0 * hx_half._duty_fixed_effectiveness(state_half),
        rel_tol=1e-9,
    )


def test_effectiveness_mode_report_metrics_includes_power_and_pressure_ratios():
    hx, state, Q_expected = _balanced_effectiveness_case(PR_hot=0.9, PR_cold=0.85)
    metrics = hx.report_metrics(state)
    assert math.isclose(metrics["power [W]"], Q_expected, rel_tol=1e-9)
    assert metrics["PR_hot [-]"] == 0.9
    assert metrics["PR_cold [-]"] == 0.85


def test_effectiveness_mode_rejects_effectiveness_out_of_range():
    with pytest.raises(ValueError, match="effectiveness"):
        HeatExchanger(name="hx1", PR_hot=1.0, PR_cold=1.0, effectiveness=1.5)


def test_effectiveness_mode_rejects_non_positive_pr():
    with pytest.raises(ValueError, match="PR_hot"):
        HeatExchanger(name="hx1", PR_hot=0.0, PR_cold=1.0, effectiveness=0.8)


# --- UA mode cross-check against MultiPassHeatExchanger --------------------


def test_ua_mode_matches_multi_pass_heat_exchanger_for_equivalent_params():
    kwargs = dict(PR_hot=0.98, PR_cold=0.97, arrangement="counterflow")
    hx = HeatExchanger(name="hx1", UA=500.0, **kwargs)
    mp = MultiPassHeatExchanger(name="hx1", UA=500.0, **kwargs)

    mdot_hot, mdot_cold = 1.0, 1.2
    P_hot, T_hot_in = 300000.0, 500.0
    P_cold, T_cold_in = 300000.0, 300.0
    h_hot_in = AIR.enthalpy_pt(P_hot, T_hot_in)
    h_cold_in = AIR.enthalpy_pt(P_cold, T_cold_in)
    state = _FakeState(
        fluid=AIR,
        mdot={
            "hx1.hot_in": mdot_hot, "hx1.hot_out": mdot_hot,
            "hx1.cold_in": mdot_cold, "hx1.cold_out": mdot_cold,
        },
        node_values={
            "hx1.hot_in": (P_hot, h_hot_in), "hx1.hot_out": (P_hot, h_hot_in),
            "hx1.cold_in": (P_cold, h_cold_in), "hx1.cold_out": (P_cold, h_cold_in),
        },
    )
    assert hx.residuals(state) == mp.residuals(state)


def test_multi_pass_heat_exchanger_is_a_heat_exchanger_subclass():
    mp = MultiPassHeatExchanger(name="hx1", UA=500.0, PR_hot=0.98, PR_cold=0.97)
    assert isinstance(mp, HeatExchanger)
    assert mp._mode == "ua"


# --- UA=None: free Newton unknown, closed by a Setpoint --------------------


def test_ua_none_is_a_free_parameter_seeded_from_ua_guess():
    hx = HeatExchanger(name="hx1", PR_hot=1.0, PR_cold=1.0, UA_guess=2.5e4)
    assert hx.free_parameters() == {"UA": 2.5e4}


def test_ua_none_end_to_end_solve_hits_setpoint_target_hot_out_temperature():
    from thermowave.components.setpoint import Setpoint
    from thermowave.components.sink import Sink
    from thermowave.components.source import Source
    from thermowave.core.network import Network

    hot_src = Source(name="hot_src", P=300000.0, T=500.0, mdot=1.0)
    cold_src = Source(name="cold_src", P=300000.0, T=300.0, mdot=1.2)
    hx = HeatExchanger(name="hx1", PR_hot=1.0, PR_cold=1.0)  # UA left free
    target_T_hot_out = 420.0
    sp = Setpoint(
        name="sp1", component=hx, free_param="UA",
        target_metric="T_hot_out [K]", value=target_T_hot_out,
    )
    hot_snk = Sink(name="hot_snk")
    cold_snk = Sink(name="cold_snk")

    network = Network(fluid=AIR)
    for c in (hot_src, cold_src, hx, sp, hot_snk, cold_snk):
        network.add_component(c)
    network.connect(hot_src, "out", hx, "hot_in")
    network.connect(hx, "hot_out", hot_snk, "in")
    network.connect(cold_src, "out", hx, "cold_in")
    network.connect(hx, "cold_out", cold_snk, "in")

    result = network.solve(tol=1e-9, max_iter=200)
    metrics = hx.report_metrics(result.state())
    assert math.isclose(metrics["T_hot_out [K]"], target_T_hot_out, abs_tol=1e-3)
    assert result.params["hx1.UA"] > 0.0
