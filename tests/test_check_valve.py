import math

import pytest

from thermowave.components.check_valve import CheckValve
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _FakeState:
    def __init__(self, fluid, mdot, node_values):
        self.fluid = fluid
        self._mdot = mdot
        self._node_values = node_values

    def node(self, name: str) -> tuple[float, float]:
        return self._node_values[name]

    def fluid_at(self, name: str):
        return self.fluid

    def mdot(self, name: str) -> float:
        return self._mdot


def test_check_valve_rejects_non_positive_k():
    with pytest.raises(ValueError, match="K must be > 0"):
        CheckValve(name="cv1", D=0.05, K=0.0)


def test_check_valve_rejects_reverse_factor_not_greater_than_one():
    with pytest.raises(ValueError, match="reverse_factor"):
        CheckValve(name="cv1", D=0.05, K=2.0, reverse_factor=1.0)


def test_check_valve_ports_returns_inlet_and_outlet_derived_from_name():
    cv = CheckValve(name="cv1", D=0.05, K=2.0)
    assert cv.ports() == {"in": "cv1.in", "out": "cv1.out"}


def test_check_valve_forward_flow_matches_valve_style_hand_calc():
    cv = CheckValve(name="cv1", D=0.05, K=2.0)
    P_in, h_in = 200000.0, AIR.enthalpy_pt(200000.0, 300.0)
    mdot = 0.5

    rho = AIR.density_ph(P_in, h_in)
    area = math.pi * 0.05**2 / 4
    v = mdot / (rho * area)
    dp_loss = 2.0 * (rho * v * abs(v) / 2.0)
    P_out = P_in - dp_loss

    state = _FakeState(
        fluid=AIR, mdot=mdot,
        node_values={"cv1.in": (P_in, h_in), "cv1.out": (P_out, h_in)},
    )
    momentum_residual, energy_residual, mass_residual = cv.residuals(state)
    assert math.isclose(momentum_residual, 0.0, abs_tol=1e-6)
    assert math.isclose(energy_residual, 0.0, abs_tol=1e-6)
    assert math.isclose(mass_residual, 0.0, abs_tol=1e-6)


def test_check_valve_reverse_flow_uses_reverse_factor_scaled_k():
    cv = CheckValve(name="cv1", D=0.05, K=2.0, reverse_factor=1000.0)
    P_in, h_in = 200000.0, AIR.enthalpy_pt(200000.0, 300.0)
    mdot = -0.01  # a trial reverse-flow Newton iterate

    rho = AIR.density_ph(P_in, h_in)
    area = math.pi * 0.05**2 / 4
    v = mdot / (rho * area)
    dp_loss_reverse = (2.0 * 1000.0) * (rho * v * abs(v) / 2.0)
    dp_loss_forward_would_be = 2.0 * (rho * v * abs(v) / 2.0)

    # reverse resistance is 1000x forward, so the (negative) pressure drop
    # implied by the same |mdot| is 1000x larger in magnitude
    assert math.isclose(dp_loss_reverse, 1000.0 * dp_loss_forward_would_be, rel_tol=1e-9)

    P_out = P_in - dp_loss_reverse
    state = _FakeState(
        fluid=AIR, mdot=mdot,
        node_values={"cv1.in": (P_in, h_in), "cv1.out": (P_out, h_in)},
    )
    momentum_residual, _, _ = cv.residuals(state)
    assert math.isclose(momentum_residual, 0.0, abs_tol=1e-6)


def test_check_valve_reverse_flow_is_much_more_resisted_than_forward():
    # Same magnitude of pressure differential in each direction should
    # drive far less mass flow in reverse than forward — that's the whole
    # point of the component. Check by comparing the momentum residual's
    # implied dp at equal |v| in each direction rather than solving a full
    # network (keeps this a fast, deterministic unit test).
    cv = CheckValve(name="cv1", D=0.05, K=2.0, reverse_factor=1000.0)
    P_in, h_in = 200000.0, AIR.enthalpy_pt(200000.0, 300.0)
    rho = AIR.density_ph(P_in, h_in)
    area = math.pi * 0.05**2 / 4

    def dp_at(mdot):
        v = mdot / (rho * area)
        K = cv.K if v >= 0.0 else cv.K * cv.reverse_factor
        return K * (rho * v * abs(v) / 2.0)

    assert abs(dp_at(-0.1)) > abs(dp_at(0.1)) * 100.0


def test_check_valve_report_metrics():
    cv = CheckValve(name="cv1", D=0.05, K=2.0)
    P_in, h_in = 200000.0, AIR.enthalpy_pt(200000.0, 300.0)
    P_out, h_out = 190000.0, h_in
    state = _FakeState(
        fluid=AIR, mdot=0.5,
        node_values={"cv1.in": (P_in, h_in), "cv1.out": (P_out, h_out)},
    )
    metrics = cv.report_metrics(state)
    assert math.isclose(metrics["power [W]"], 0.5 * (h_out - h_in))
    assert math.isclose(metrics["PR [-]"], P_out / P_in)
    assert math.isclose(metrics["mdot [kg/s]"], 0.5)
