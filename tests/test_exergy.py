import math

import pytest

from thermowave.components.compressor import Compressor
from thermowave.components.electric_motor import ElectricMotor
from thermowave.components.multi_pass_heat_exchanger import MultiPassHeatExchanger
from thermowave.components.pump import Pump
from thermowave.components.simple_compressor import SimpleCompressor
from thermowave.components.simple_condenser import SimpleCondenser
from thermowave.components.simple_heat_exchanger import SimpleHeatExchanger
from thermowave.components.simple_turbine import SimpleTurbine
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.exergy import (
    ComponentExergyCost,
    exergy_report,
    node_exergy,
)
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)
T0, P0 = 298.15, 101325.0


class _FakeState:
    def __init__(self, fluid, node_values, mdots=None, params=None):
        self.fluid = fluid
        self._node_values = node_values
        self._mdots = mdots or {}
        self._params = params or {}

    def node(self, name):
        return self._node_values[name]

    def fluid_at(self, name):
        return self.fluid

    def mdot(self, name):
        if name not in self._mdots:
            raise KeyError(name)
        return self._mdots[name]

    def param(self, name):
        return self._params[name]


class _FakeResult:
    """Minimal SolveResult-shaped stub for exergy_report(): only the fields
    it actually reads (components, node_order, .state())."""

    def __init__(self, components, node_order, state):
        self.components = components
        self.node_order = node_order
        self._state = state

    def state(self):
        return self._state


# --- node_exergy -------------------------------------------------------


def test_node_exergy_is_zero_at_the_dead_state():
    h0 = AIR.enthalpy_pt(P0, T0)
    state = _FakeState(fluid=AIR, node_values={"n": (P0, h0)}, mdots={"n": 1.0})
    ex = node_exergy(state, "n", T0, P0)
    assert math.isclose(ex.e_PH, 0.0, abs_tol=1e-6)
    assert math.isclose(ex.e_T, 0.0, abs_tol=1e-6)
    assert math.isclose(ex.e_M, 0.0, abs_tol=1e-6)


def test_node_exergy_e_t_plus_e_m_equals_e_ph():
    P, T = 300000.0, 450.0
    h = AIR.enthalpy_pt(P, T)
    state = _FakeState(fluid=AIR, node_values={"n": (P, h)}, mdots={"n": 2.0})
    ex = node_exergy(state, "n", T0, P0)
    assert math.isclose(ex.e_T + ex.e_M, ex.e_PH, rel_tol=1e-9)


def test_node_exergy_e_ph_matches_direct_two_state_formula():
    P, T = 300000.0, 450.0
    h = AIR.enthalpy_pt(P, T)
    h0 = AIR.enthalpy_pt(P0, T0)
    s = AIR.entropy_ph(P, h)
    s0 = AIR.entropy_ph(P0, h0)
    expected = (h - h0) - T0 * (s - s0)

    state = _FakeState(fluid=AIR, node_values={"n": (P, h)}, mdots={"n": 2.0})
    ex = node_exergy(state, "n", T0, P0)
    assert math.isclose(ex.e_PH, expected, rel_tol=1e-9)


def test_node_exergy_flow_rate_scales_by_mdot():
    P, T = 300000.0, 450.0
    h = AIR.enthalpy_pt(P, T)
    state = _FakeState(fluid=AIR, node_values={"n": (P, h)}, mdots={"n": 3.0})
    ex = node_exergy(state, "n", T0, P0)
    assert math.isclose(ex.E_PH, 3.0 * ex.e_PH, rel_tol=1e-9)


def test_node_exergy_flow_rate_is_none_without_mdot():
    P, T = 300000.0, 450.0
    h = AIR.enthalpy_pt(P, T)
    state = _FakeState(fluid=AIR, node_values={"n": (P, h)})  # no mdots given
    ex = node_exergy(state, "n", T0, P0)
    assert ex.E_PH is None
    assert ex.E_T is None
    assert ex.E_M is None


def test_node_exergy_raises_for_cantera_fluid():
    ct = pytest.importorskip("cantera")
    from thermowave.fluids.cantera_fluid import CanteraFluid

    gas = CanteraFluid(name="air", composition={"N2": 0.79, "O2": 0.21})
    P, T = 101325.0, 400.0
    h = gas.enthalpy_pt(P, T)
    state = _FakeState(fluid=gas, node_values={"n": (P, h)}, mdots={"n": 1.0})
    with pytest.raises(ValueError, match="entropy_ph"):
        node_exergy(state, "n", T0, P0)


# --- per-component-type costing (via a real Network solve for realistic states) --


def _solved_state(components, connections, fluid=AIR):
    net = Network(fluid=fluid)
    for c in components:
        net.add_component(c)
    for a, a_port, b, b_port in connections:
        net.connect(a, a_port, b, b_port)
    result = net.solve(progress=False)
    assert result.converged
    return result


def test_expander_cost_matches_hand_calc():
    src = Source(name="src", P=500000.0, T=600.0, mdot=2.0)
    turb = SimpleTurbine(name="turb", PR=4.0, eta_s=0.85, gamma=1.4)
    snk = Sink(name="snk")
    result = _solved_state([src, turb, snk], [(src, "out", turb, "in"), (turb, "out", snk, "in")])

    state = result.state()
    from thermowave.core.exergy import _cost_component

    cost = _cost_component(turb, state, frozenset(), {}, T0, P0)
    E_in = node_exergy(state, "turb.in", T0, P0).E_PH
    E_out = node_exergy(state, "turb.out", T0, P0).E_PH
    expected_E_F = E_in - E_out
    expected_E_P = turb.report_metrics(state)["power [W]"]
    assert isinstance(cost, ComponentExergyCost)
    assert math.isclose(cost.E_F, expected_E_F, rel_tol=1e-9)
    assert math.isclose(cost.E_P, expected_E_P, rel_tol=1e-9)
    assert math.isclose(cost.E_D, expected_E_F - expected_E_P, rel_tol=1e-9)
    assert cost.E_D >= -1e-6  # irreversible (eta_s < 1) expander destroys exergy


def test_compressor_cost_has_swapped_fuel_product_roles():
    src = Source(name="src", P=101325.0, T=300.0, mdot=2.0)
    comp = SimpleCompressor(name="comp", PR=4.0, eta_s=0.85, gamma=1.4)
    snk = Sink(name="snk")
    result = _solved_state([src, comp, snk], [(src, "out", comp, "in"), (comp, "out", snk, "in")])

    state = result.state()
    from thermowave.core.exergy import _cost_component

    cost = _cost_component(comp, state, frozenset(), {}, T0, P0)
    expected_E_F = comp.report_metrics(state)["power [W]"]
    expected_E_P = (
        node_exergy(state, "comp.out", T0, P0).E_PH - node_exergy(state, "comp.in", T0, P0).E_PH
    )
    assert math.isclose(cost.E_F, expected_E_F, rel_tol=1e-9)
    assert math.isclose(cost.E_P, expected_E_P, rel_tol=1e-9)
    assert cost.E_D >= -1e-6


def test_two_stream_hx_default_costs_cold_side_gain_as_product():
    hot_src = Source(name="hot_src", P=101325.0, T=500.0, mdot=1.0)
    cold_src = Source(name="cold_src", P=101325.0, T=300.0, mdot=1.0)
    hx = SimpleHeatExchanger(name="hx", effectiveness=0.7, PR_hot=1.0, PR_cold=1.0)
    hot_snk = Sink(name="hot_snk")
    cold_snk = Sink(name="cold_snk")
    result = _solved_state(
        [hot_src, cold_src, hx, hot_snk, cold_snk],
        [
            (hot_src, "out", hx, "hot_in"), (hx, "hot_out", hot_snk, "in"),
            (cold_src, "out", hx, "cold_in"), (hx, "cold_out", cold_snk, "in"),
        ],
    )
    state = result.state()
    from thermowave.core.exergy import _cost_component

    cost = _cost_component(hx, state, frozenset(), {}, T0, P0)
    expected_E_F = (
        node_exergy(state, "hx.hot_in", T0, P0).E_PH - node_exergy(state, "hx.hot_out", T0, P0).E_PH
    )
    expected_E_P = (
        node_exergy(state, "hx.cold_out", T0, P0).E_PH - node_exergy(state, "hx.cold_in", T0, P0).E_PH
    )
    assert math.isclose(cost.E_F, expected_E_F, rel_tol=1e-9)
    assert math.isclose(cost.E_P, expected_E_P, rel_tol=1e-9)


def test_two_stream_hx_dissipative_has_no_product_and_feeds_e_l_not_e_d():
    hot_src = Source(name="hot_src", P=101325.0, T=500.0, mdot=1.0)
    cold_src = Source(name="cold_src", P=101325.0, T=300.0, mdot=1.0)
    hx = SimpleHeatExchanger(name="cooler", effectiveness=0.7, PR_hot=1.0, PR_cold=1.0)
    hot_snk = Sink(name="hot_snk")
    cold_snk = Sink(name="cold_snk")
    net = Network(fluid=AIR)
    for c in (hot_src, cold_src, hx, hot_snk, cold_snk):
        net.add_component(c)
    net.connect(hot_src, "out", hx, "hot_in")
    net.connect(hx, "hot_out", hot_snk, "in")
    net.connect(cold_src, "out", hx, "cold_in")
    net.connect(hx, "cold_out", cold_snk, "in")
    result = net.solve(progress=False)
    assert result.converged

    report = exergy_report(result, T0, P0, fuel=[0.0], product=[0.0], dissipative=["cooler"])
    cost = report.component_costs["cooler"]
    assert cost.E_P is None
    assert cost.E_D is None
    assert report.E_L > 0.0
    assert math.isclose(report.E_L, cost.E_F, rel_tol=1e-9)
    assert report.E_D == 0.0  # nothing else in this tiny network is costed


def test_dissipative_rejects_non_two_stream_hx_name():
    src = Source(name="src", P=101325.0, T=400.0, mdot=1.0)
    snk = Sink(name="snk")
    net = Network(fluid=AIR)
    net.add_component(src)
    net.add_component(snk)
    net.connect(src, "out", snk, "in")
    result = net.solve(progress=False)
    with pytest.raises(ValueError, match="SimpleHeatExchanger"):
        exergy_report(result, T0, P0, fuel=[0.0], product=[0.0], dissipative=["src"])


def test_single_stream_condenser_costed_only_with_source_temperature():
    pytest.importorskip("CoolProp")
    from thermowave.fluids.real_fluid import CoolPropFluid

    water = CoolPropFluid(name="Water")
    # Operate above ambient pressure (P0=101325) so the hybrid (T, P0) state
    # node_exergy evaluates internally doesn't land exactly on water's own
    # saturation dome at 1 atm (373.12 K) -- a real CoolProp edge case, not
    # an exergy-formula bug: T at this condenser's own operating pressure
    # is comfortably above that coincidence.
    P_operating = 200000.0
    T_sat = water.saturation_temperature(P_operating)
    src = Source(name="src", P=P_operating, T=T_sat + 0.01, mdot=1.0)
    cond = SimpleCondenser(name="cond", outlet_quality=0.0)
    snk = Sink(name="snk")
    result = _solved_state(
        [src, cond, snk], [(src, "out", cond, "in"), (cond, "out", snk, "in")], fluid=water
    )
    state = result.state()
    from thermowave.core.exergy import _cost_component

    assert _cost_component(cond, state, frozenset(), {}, T0, P0) is None

    T_source = T_sat - 20.0  # a cooling-water sink, below the condenser saturation temp
    cost = _cost_component(cond, state, frozenset(), {"cond": T_source}, T0, P0)
    Q = cond.report_metrics(state)["power [W]"]
    expected_E_F = abs(Q) * (1.0 - T0 / T_source)
    assert math.isclose(cost.E_F, expected_E_F, rel_tol=1e-9)


def test_power_only_components_are_not_costed():
    motor = ElectricMotor(name="motor", component=SimpleCompressor(
        name="c", PR=2.0, eta_s=0.8, gamma=1.4
    ), efficiency=0.95)
    from thermowave.core.exergy import _cost_component

    class _FakeStateForMotor:
        pass

    cost = _cost_component(motor, _FakeStateForMotor(), frozenset(), {}, T0, P0)
    assert cost is None


# --- network-level totals ------------------------------------------------


def test_network_totals_use_callable_fuel_and_product_like_setpoint():
    src = Source(name="src", P=500000.0, T=600.0, mdot=2.0)
    turb = SimpleTurbine(name="turb", PR=4.0, eta_s=0.85, gamma=1.4)
    snk = Sink(name="snk")
    result = _solved_state([src, turb, snk], [(src, "out", turb, "in"), (turb, "out", snk, "in")])

    report = exergy_report(
        result, T0, P0,
        fuel=[lambda s: node_exergy(s, "src.out", T0, P0).E_PH],
        product=[lambda s: turb.report_metrics(s)["power [W]"]],
    )
    expected_fuel = node_exergy(result.state(), "src.out", T0, P0).E_PH
    expected_product = turb.report_metrics(result.state())["power [W]"]
    assert math.isclose(report.E_F, expected_fuel, rel_tol=1e-9)
    assert math.isclose(report.E_P, expected_product, rel_tol=1e-9)
    assert math.isclose(report.epsilon, expected_product / expected_fuel, rel_tol=1e-9)


# --- end-to-end energy-conservation sanity check -------------------------


def test_end_to_end_turbine_and_compressor_satisfy_e_f_minus_e_p_equals_e_d():
    src = Source(name="src", P=101325.0, T=300.0, mdot=2.0)
    comp = SimpleCompressor(name="comp", PR=4.0, eta_s=0.85, gamma=1.4)
    turb = SimpleTurbine(name="turb", PR=3.5, eta_s=0.88, gamma=1.4)
    snk = Sink(name="snk")
    result = _solved_state(
        [src, comp, turb, snk],
        [(src, "out", comp, "in"), (comp, "out", turb, "in"), (turb, "out", snk, "in")],
    )

    report = exergy_report(result, T0, P0, fuel=[0.0], product=[0.0])
    for name in ("comp", "turb"):
        cost = report.component_costs[name]
        assert math.isclose(cost.E_F - cost.E_P, cost.E_D, abs_tol=1e-6)
        assert cost.E_D >= -1e-6  # both are irreversible (eta_s < 1)
