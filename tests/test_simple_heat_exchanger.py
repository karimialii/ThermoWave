import math

from thermowave.components.simple_heat_exchanger import SimpleHeatExchanger
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _FakeState:
    def __init__(self, fluid, mdot, node_values, params=None):
        self.fluid = fluid
        self._mdot = mdot
        self._node_values = node_values
        self._params = params or {}

    def node(self, name: str) -> tuple[float, float]:
        return self._node_values[name]

    def fluid_at(self, name: str):
        return self.fluid

    def mdot(self, name: str) -> float:
        return self._mdot[name]

    def param(self, name: str) -> float:
        return self._params[name]


def _case(Q=50_000.0, PR=0.98):
    mdot = 1.5
    P_in, T_in = 300000.0, 400.0
    h_in = AIR.enthalpy_pt(P_in, T_in)

    hx = SimpleHeatExchanger(name="hx1", Q=Q, PR=PR)

    h_out = h_in + Q / mdot
    P_out = PR * P_in

    state = _FakeState(
        fluid=AIR,
        mdot={"hx1.in": mdot, "hx1.out": mdot},
        node_values={"hx1.in": (P_in, h_in), "hx1.out": (P_out, h_out)},
    )
    return hx, state, Q


def test_ports_returns_in_and_out_derived_from_name():
    hx = SimpleHeatExchanger(name="hx1", Q=1000.0)
    assert hx.ports() == {"in": "hx1.in", "out": "hx1.out"}


def test_residuals_all_zero_at_exact_hand_calc_solution():
    hx, state, _Q = _case()
    residuals = hx.residuals(state)
    assert len(residuals) == 3
    for r in residuals:
        assert math.isclose(r, 0.0, abs_tol=1e-6)


def test_negative_q_cools_the_stream():
    hx, state, Q = _case(Q=-40_000.0)
    assert Q < 0.0
    momentum, energy, mass = hx.residuals(state)
    assert math.isclose(energy, 0.0, abs_tol=1e-6)


def test_momentum_residual_matches_pressure_ratio_hand_calc():
    hx, state, _Q = _case(PR=0.9)
    momentum_residual, _energy, _mass = hx.residuals(state)
    assert math.isclose(momentum_residual, 0.0, abs_tol=1e-6)


def test_report_metrics_includes_power_and_pressure_ratio():
    hx, state, Q_expected = _case(PR=0.85)
    metrics = hx.report_metrics(state)
    assert math.isclose(metrics["power [W]"], Q_expected, rel_tol=1e-9)
    assert metrics["PR [-]"] == 0.85


def test_simple_heat_exchanger_rejects_non_positive_pr():
    try:
        SimpleHeatExchanger(name="hx1", Q=1000.0, PR=0.0)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_free_parameters_empty_when_q_given():
    hx = SimpleHeatExchanger(name="hx1", Q=1000.0)
    assert hx.free_parameters() == {}


def test_free_parameters_includes_q_when_q_omitted():
    hx = SimpleHeatExchanger(name="hx1", Q=None)
    assert hx.free_parameters() == {"Q": 0.0}


def test_residuals_use_q_from_state_param_when_q_omitted():
    mdot = 1.5
    P_in, T_in = 300000.0, 400.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    Q = 25_000.0
    hx = SimpleHeatExchanger(name="hx1", Q=None, PR=1.0)
    h_out = h_in + Q / mdot
    state = _FakeState(
        fluid=AIR,
        mdot={"hx1.in": mdot, "hx1.out": mdot},
        node_values={"hx1.in": (P_in, h_in), "hx1.out": (P_in, h_out)},
        params={"hx1.Q": Q},
    )
    momentum, energy, mass = hx.residuals(state)
    assert math.isclose(energy, 0.0, abs_tol=1e-6)
