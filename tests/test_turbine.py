import math
from pathlib import Path

from thermowave.components.turbine import Turbine
from thermowave.fluids.ideal_gas import IdealGasFluid

_MAP_PATH = str(Path(__file__).parent / "fixtures" / "simple_turbine_map.tur")

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _FakeHeatPath:
    def __init__(self, Q_value):
        self._Q_value = Q_value

    def Q(self, state):
        return self._Q_value


class _FakeState:
    """Minimal stand-in for core.network.NetworkState."""

    def __init__(self, fluid, mdot, node_values: dict[str, tuple[float, float]], params=None):
        self.fluid = fluid
        self._mdot = mdot
        self._node_values = node_values
        self._params = params or {}

    def node(self, name: str) -> tuple[float, float]:
        return self._node_values[name]

    def fluid_at(self, name: str):
        return self.fluid

    def mdot(self, name: str) -> float:
        return self._mdot

    def param(self, name: str) -> float:
        return self._params[name]


def test_turbine_ports_returns_inlet_and_outlet_derived_from_name():
    turb = Turbine(name="t1", map_path=_MAP_PATH, N=3000.0, gamma=1.4)
    assert turb.ports() == {"in": "t1.in", "out": "t1.out"}


def test_turbine_residuals_matches_map_lookup_hand_calc():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1.4
    P_in, T_in = 400000.0, 1000.0
    h_in = air.enthalpy_pt(P_in, T_in)

    # Corrected speed A = 50 [rev/s/sqrt(K)] exactly, at T_in = 1000 K.
    N = 50.0 * math.sqrt(T_in) * 60.0
    turb = Turbine(name="t1", map_path=_MAP_PATH, N=N, gamma=gamma)

    # Corrected mass flow B = 2 exactly, at P_in = 400000 Pa = 4 Bar.
    mdot = 2.0 * (P_in / 1.0e5) / math.sqrt(T_in)

    PR, eta_s = 2.0, 0.86  # map values at A=50, B=2 (see fixture)
    P_out_guess = P_in / PR
    T_out_isentropic = T_in * (1.0 / PR) ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_out_guess, T_out_isentropic)
    dh_actual = eta_s * (h_in - h_out_isentropic)
    h_out_guess = h_in - dh_actual

    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={
            "t1.in": (P_in, h_in),
            "t1.out": (P_out_guess, h_out_guess),
        },
    )
    momentum_residual, energy_residual, mass_residual = turb.residuals(state)
    assert math.isclose(momentum_residual, 0.0, abs_tol=1e-6)
    assert math.isclose(energy_residual, 0.0, abs_tol=1e-3)
    assert math.isclose(mass_residual, 0.0, abs_tol=1e-6)


def test_turbine_gamma_none_derives_from_fluid_and_matches_explicit_gamma():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1005.0 / (1005.0 - 287.05)
    P_in, T_in = 400000.0, 1000.0
    h_in = air.enthalpy_pt(P_in, T_in)
    N = 50.0 * math.sqrt(T_in) * 60.0
    mdot = 2.0 * (P_in / 1.0e5) / math.sqrt(T_in)
    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={
            "t1.in": (P_in, h_in),
            "t1.out": (P_in / 2.0, h_in),  # arbitrary guess, same for both
        },
    )

    explicit = Turbine(name="t1", map_path=_MAP_PATH, N=N, gamma=gamma)
    derived = Turbine(name="t1", map_path=_MAP_PATH, N=N, gamma=None)
    assert derived.gamma is None

    assert derived.residuals(state) == explicit.residuals(state)


def test_turbine_factor_overrides_defaults_to_file_values():
    turb_default = Turbine(name="t1", map_path=_MAP_PATH, N=3000.0, gamma=1.4)
    turb_explicit_none = Turbine(
        name="t1", map_path=_MAP_PATH, N=3000.0, gamma=1.4, factor_overrides=None
    )
    A, B = turb_default.map.mid_speed(), 1.0
    assert math.isclose(
        turb_default.map.pressure_ratio(A, B), turb_explicit_none.map.pressure_ratio(A, B)
    )


def test_turbine_factor_overrides_calibrates_pressure_ratio():
    turb_default = Turbine(name="t1", map_path=_MAP_PATH, N=3000.0, gamma=1.4)
    turb_calibrated = Turbine(
        name="t1",
        map_path=_MAP_PATH,
        N=3000.0,
        gamma=1.4,
        factor_overrides={"C_fact": 1.05},
    )
    A, B = turb_default.map.mid_speed(), 1.0
    assert math.isclose(
        turb_calibrated.map.pressure_ratio(A, B),
        1.05 * turb_default.map.pressure_ratio(A, B),
        rel_tol=1e-9,
    )


def test_turbine_heat_path_defaults_to_none():
    turb = Turbine(name="t1", map_path=_MAP_PATH, N=3000.0, gamma=1.4)
    assert turb.heat_path is None


def test_turbine_heat_path_adds_q_over_mdot_to_energy_residual():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1.4
    P_in, T_in = 400000.0, 1000.0
    h_in = air.enthalpy_pt(P_in, T_in)
    N = 50.0 * math.sqrt(T_in) * 60.0
    mdot = 2.0 * (P_in / 1.0e5) / math.sqrt(T_in)

    PR, eta_s = 2.0, 0.86  # map values at A=50, B=2 (see fixture)
    P_out = P_in / PR
    T_out_isentropic = T_in * (1.0 / PR) ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_out, T_out_isentropic)
    dh_actual = eta_s * (h_in - h_out_isentropic)
    # Same h_out guess the no-heat-loss case converges to (residual ~ 0 there).
    h_out_guess = h_in - dh_actual

    state = _FakeState(
        fluid=air, mdot=mdot, node_values={"t1.in": (P_in, h_in), "t1.out": (P_out, h_out_guess)},
    )

    baseline = Turbine(name="t1", map_path=_MAP_PATH, N=N, gamma=gamma)
    _, baseline_energy_residual, _ = baseline.residuals(state)

    Q = 5000.0
    lossy = Turbine(name="t1", map_path=_MAP_PATH, N=N, gamma=gamma)
    lossy.heat_path = _FakeHeatPath(Q_value=Q)
    _, lossy_energy_residual, _ = lossy.residuals(state)

    assert math.isclose(
        lossy_energy_residual - baseline_energy_residual, Q / mdot, rel_tol=1e-9
    )


def test_turbine_report_metrics_includes_q_loss():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    P_in, T_in = 400000.0, 1000.0
    h_in = air.enthalpy_pt(P_in, T_in)
    N = 50.0 * math.sqrt(T_in) * 60.0
    mdot = 2.0 * (P_in / 1.0e5) / math.sqrt(T_in)
    state = _FakeState(
        fluid=air, mdot=mdot, node_values={"t1.in": (P_in, h_in), "t1.out": (P_in / 2.0, h_in)},
    )

    turb = Turbine(name="t1", map_path=_MAP_PATH, N=N, gamma=1.4)
    assert turb.report_metrics(state)["Q_loss [W]"] == 0.0

    turb.heat_path = _FakeHeatPath(Q_value=1234.0)
    assert turb.report_metrics(state)["Q_loss [W]"] == 1234.0


def test_turbine_report_metrics_reflects_map_lookup():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1.4
    P_in, T_in = 400000.0, 1000.0
    h_in = air.enthalpy_pt(P_in, T_in)

    N = 70.0 * math.sqrt(T_in) * 60.0
    turb = Turbine(name="t1", map_path=_MAP_PATH, N=N, gamma=gamma)
    mdot = 2.0 * (P_in / 1.0e5) / math.sqrt(T_in)

    PR, eta_s = 2.4, 0.88  # map values at A=70, B=2 (see fixture)
    P_out = P_in / PR
    T_out_isentropic = T_in * (1.0 / PR) ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_out, T_out_isentropic)
    h_out = h_in - eta_s * (h_in - h_out_isentropic)

    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={
            "t1.in": (P_in, h_in),
            "t1.out": (P_out, h_out),
        },
    )
    metrics = turb.report_metrics(state)
    assert math.isclose(metrics["PR [-]"], PR, rel_tol=1e-9)
    assert math.isclose(metrics["eta_s [-]"], eta_s, rel_tol=1e-9)
    assert math.isclose(metrics["power [W]"], mdot * (h_in - h_out), rel_tol=1e-9)


def test_turbine_free_parameters_empty_when_n_given():
    turb = Turbine(name="t1", map_path=_MAP_PATH, N=3000.0, gamma=1.4)
    assert turb.free_parameters() == {}


def test_turbine_free_parameters_includes_n_when_n_omitted():
    turb = Turbine(name="t1", map_path=_MAP_PATH, gamma=1.4)
    params = turb.free_parameters()
    assert set(params.keys()) == {"N"}
    assert params["N"] > 0


def test_turbine_residuals_uses_n_from_state_param_when_n_omitted():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1.4
    P_in, T_in = 400000.0, 1000.0
    h_in = air.enthalpy_pt(P_in, T_in)

    N = 50.0 * math.sqrt(T_in) * 60.0
    mdot = 2.0 * (P_in / 1.0e5) / math.sqrt(T_in)

    turb = Turbine(name="t1", map_path=_MAP_PATH, gamma=gamma)  # N left free
    PR, eta_s = 2.0, 0.86  # map values at A=50, B=2 (see fixture)
    P_out = P_in / PR
    T_out_isentropic = T_in * (1.0 / PR) ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_out, T_out_isentropic)
    h_out = h_in - eta_s * (h_in - h_out_isentropic)

    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={
            "t1.in": (P_in, h_in),
            "t1.out": (P_out, h_out),
        },
        params={"t1.N": N},
    )
    momentum_residual, energy_residual, mass_residual = turb.residuals(state)
    assert math.isclose(momentum_residual, 0.0, abs_tol=1e-6)
    assert math.isclose(energy_residual, 0.0, abs_tol=1e-3)
    assert math.isclose(mass_residual, 0.0, abs_tol=1e-6)


def test_turbine_rejects_gamma_not_greater_than_one():
    try:
        Turbine(name="t1", map_path=_MAP_PATH, gamma=0.9)
        assert False, "expected ValueError"
    except ValueError:
        pass
