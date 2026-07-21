import math

import pytest

from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.tank import Tank
from thermowave.components.valve import Valve
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _FakeState:
    """Minimal stand-in for core.network.NetworkState."""

    def __init__(self, fluid, node_values, mdots, params):
        self.fluid = fluid
        self._node_values = node_values
        self._mdots = mdots
        self._params = params

    def node(self, name: str) -> tuple[float, float]:
        return self._node_values[name]

    def fluid_at(self, name: str):
        return self.fluid

    def mdot(self, name: str) -> float:
        return self._mdots[name]

    def param(self, name: str) -> float:
        return self._params[name]


def test_tank_rejects_non_positive_volume():
    with pytest.raises(ValueError, match="V must be > 0"):
        Tank(name="t1", V=0.0, P0=101325.0, T0=300.0, fluid=AIR)


def test_tank_ports_returns_inlet_and_outlet_derived_from_name():
    tank = Tank(name="t1", V=0.5, P0=101325.0, T0=300.0, fluid=AIR)
    assert tank.ports() == {"in": "t1.in", "out": "t1.out"}


def test_tank_differential_parameters_seeded_from_p0_t0():
    tank = Tank(name="t1", V=0.5, P0=101325.0, T0=300.0, fluid=AIR)
    params = tank.differential_parameters()
    assert math.isclose(params["P"], 101325.0)
    assert math.isclose(params["h"], AIR.enthalpy_pt(101325.0, 300.0))


def test_tank_residuals_ties_outlet_node_to_tank_state():
    tank = Tank(name="t1", V=0.5, P0=101325.0, T0=300.0, fluid=AIR)
    h0 = AIR.enthalpy_pt(101325.0, 300.0)
    state = _FakeState(
        fluid=AIR,
        node_values={"t1.in": (150000.0, 320000.0), "t1.out": (140000.0, 300000.0)},
        mdots={"t1.in": 0.5, "t1.out": 0.5},
        params={"t1.P": 101325.0, "t1.h": h0},
    )
    momentum_residual, energy_residual = tank.residuals(state)
    assert math.isclose(momentum_residual, 140000.0 - 101325.0)
    assert math.isclose(energy_residual, 300000.0 - h0)


def test_tank_state_derivative_is_zero_at_matched_inflow_outflow_and_no_heat_loss():
    tank = Tank(name="t1", V=0.5, P0=150000.0, T0=300.0, fluid=AIR)
    h0 = AIR.enthalpy_pt(150000.0, 300.0)
    state = _FakeState(
        fluid=AIR,
        node_values={"t1.in": (150000.0, h0), "t1.out": (150000.0, h0)},
        mdots={"t1.in": 0.5, "t1.out": 0.5},
        params={"t1.P": 150000.0, "t1.h": h0},
    )
    derivative = tank.state_derivative(state)
    assert math.isclose(derivative["P"], 0.0, abs_tol=1e-6)
    assert math.isclose(derivative["h"], 0.0, abs_tol=1e-6)


def test_tank_pressurizes_when_filling_with_no_outflow():
    tank = Tank(name="t1", V=0.5, P0=150000.0, T0=300.0, fluid=AIR)
    h0 = AIR.enthalpy_pt(150000.0, 300.0)
    h_in = AIR.enthalpy_pt(200000.0, 320.0)
    state = _FakeState(
        fluid=AIR,
        node_values={"t1.in": (200000.0, h_in), "t1.out": (150000.0, h0)},
        mdots={"t1.in": 0.3, "t1.out": 0.0},
        params={"t1.P": 150000.0, "t1.h": h0},
    )
    derivative = tank.state_derivative(state)
    assert derivative["P"] > 0.0  # mass flowing in with nowhere to go -> pressure rises


def test_tank_state_derivative_matches_hand_solved_2x2_system():
    tank = Tank(name="t1", V=0.5, P0=150000.0, T0=300.0, fluid=AIR)
    P, h = 150000.0, AIR.enthalpy_pt(150000.0, 300.0)
    h_in = AIR.enthalpy_pt(200000.0, 320.0)
    mdot_in, mdot_out = 0.3, 0.1
    state = _FakeState(
        fluid=AIR,
        node_values={"t1.in": (200000.0, h_in), "t1.out": (P, h)},
        mdots={"t1.in": mdot_in, "t1.out": mdot_out},
        params={"t1.P": P, "t1.h": h},
    )

    rho = AIR.density_ph(P, h)
    eps_P = max(abs(P) * 1.0e-6, 1.0)
    eps_h = max(abs(h) * 1.0e-6, 1.0)
    drho_dP = (AIR.density_ph(P + eps_P, h) - rho) / eps_P
    drho_dh = (AIR.density_ph(P, h + eps_h) - rho) / eps_h

    a11, a12, b1 = 0.5 * drho_dP, 0.5 * drho_dh, mdot_in - mdot_out
    a21, a22, b2 = -0.5, rho * 0.5, mdot_in * (h_in - h)
    det = a11 * a22 - a12 * a21
    expected_dPdt = (b1 * a22 - a12 * b2) / det
    expected_dhdt = (a11 * b2 - a21 * b1) / det

    derivative = tank.state_derivative(state)
    assert math.isclose(derivative["P"], expected_dPdt, rel_tol=1e-9)
    assert math.isclose(derivative["h"], expected_dhdt, rel_tol=1e-9)


def test_tank_report_metrics_reflects_state_and_mdots():
    tank = Tank(name="t1", V=0.5, P0=150000.0, T0=300.0, fluid=AIR)
    h = AIR.enthalpy_pt(150000.0, 320.0)
    state = _FakeState(
        fluid=AIR,
        node_values={"t1.in": (150000.0, h), "t1.out": (150000.0, h)},
        mdots={"t1.in": 0.4, "t1.out": 0.3},
        params={"t1.P": 150000.0, "t1.h": h},
    )
    metrics = tank.report_metrics(state)
    assert math.isclose(metrics["P [Pa]"], 150000.0)
    assert math.isclose(metrics["T [K]"], 320.0, rel_tol=1e-6)
    assert math.isclose(metrics["mdot_in [kg/s]"], 0.4)
    assert math.isclose(metrics["mdot_out [kg/s]"], 0.3)
    assert math.isclose(metrics["V [m^3]"], 0.5)


def test_tank_steady_state_solve_balances_inflow_and_outflow():
    src = Source(name="src", P=150000.0, T=300.0, mdot=0.5)
    tank = Tank(name="tank", V=0.5, P0=150000.0, T0=300.0, fluid=AIR)
    valve = Valve(name="valve", D=0.05, K=2.0, opening=1.0)
    snk = Sink(name="snk", P=101325.0)

    network = Network(fluid=AIR)
    for component in (src, tank, valve, snk):
        network.add_component(component)
    network.connect(src, "out", tank, "in")
    network.connect(tank, "out", valve, "in")
    network.connect(valve, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=200, verbose=False)
    assert result.converged
    assert math.isclose(result.node_mdot["tank.in"], result.node_mdot["tank.out"], rel_tol=1e-6)


def test_tank_transient_fills_toward_new_steady_state_after_inflow_step():
    src = Source(name="src", P=150000.0, T=300.0, mdot=0.5)
    tank = Tank(name="tank", V=0.5, P0=101325.0, T0=300.0, fluid=AIR)
    valve = Valve(name="valve", D=0.05, K=2.0, opening=1.0)
    snk = Sink(name="snk", P=101325.0)

    network = Network(fluid=AIR)
    for component in (src, tank, valve, snk):
        network.add_component(component)
    network.connect(src, "out", tank, "in")
    network.connect(tank, "out", valve, "in")
    network.connect(valve, "out", snk, "in")

    initial = network.solve(tol=1e-8, max_iter=200)
    P_start = initial.params["tank.P"]

    src.mdot = 0.8  # step up inflow -> tank should re-pressurize over time
    history = network.solve_transient(
        duration=5.0, dt=0.05, initial=initial, tol=1e-6, max_iter=300, damping=0.5,
    )
    P_mid = history.steps[len(history.steps) // 2].params["tank.P"]
    P_end = history.steps[-1].params["tank.P"]

    # pressure rises monotonically toward the new (higher) steady state, and
    # outflow lags inflow throughout (real accumulation, not instant mixing)
    assert P_start < P_mid < P_end
    for step in history.steps[1:]:
        assert step.node_mdot["tank.out"] < step.node_mdot["tank.in"]
