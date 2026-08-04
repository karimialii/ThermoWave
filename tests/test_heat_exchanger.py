import math

import pytest

from thermowave.components.heat_exchanger import HeatExchanger, MultiPassHeatExchanger
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


# --- construction / mode selection ----------------------------------------


def test_rejects_both_effectiveness_and_ua():
    with pytest.raises(ValueError, match="exactly one"):
        HeatExchanger(name="hx1", PR_hot=1.0, PR_cold=1.0, effectiveness=0.7, UA=500.0)


def test_rejects_neither_effectiveness_nor_ua():
    with pytest.raises(ValueError, match="exactly one"):
        HeatExchanger(name="hx1", PR_hot=1.0, PR_cold=1.0)


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
