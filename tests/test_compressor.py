import math
from pathlib import Path

from thermowave.components.compressor import Compressor
from thermowave.fluids.ideal_gas import IdealGasFluid

_MAP_PATH = str(Path(__file__).parent / "fixtures" / "simple_compressor_map.cop")

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


def test_compressor_ports_returns_inlet_and_outlet_derived_from_name():
    comp = Compressor(name="c1", map_path=_MAP_PATH, N=3000.0, gamma=1.4)
    assert comp.ports() == {"in": "c1.in", "out": "c1.out"}


def test_compressor_residuals_matches_map_lookup_hand_calc():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1.4
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)
    mdot = 2.0

    # Corrected speed A = 50 [rev/s/sqrt(K)] exactly, at T_in = 300 K.
    N = 50.0 * math.sqrt(T_in) * 60.0
    comp = Compressor(name="c1", map_path=_MAP_PATH, N=N, gamma=gamma)

    # Corrected mass flow B = 2 exactly, at P_in = 101325 Pa = 1.01325 Bar.
    mdot = 2.0 * (P_in / 1.0e5) / math.sqrt(T_in)

    PR, eta_s = 3.0, 0.8  # map values at A=50, B=2 (see fixture)
    T_out_isentropic = T_in * PR ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_in * PR, T_out_isentropic)
    dh_actual = (h_out_isentropic - h_in) / eta_s
    P_out_guess = PR * P_in
    h_out_guess = h_in + dh_actual

    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={
            "c1.in": (P_in, h_in),
            "c1.out": (P_out_guess, h_out_guess),
        },
    )
    momentum_residual, energy_residual, mass_residual = comp.residuals(state)
    assert math.isclose(momentum_residual, 0.0, abs_tol=1e-6)
    assert math.isclose(energy_residual, 0.0, abs_tol=1e-3)
    assert math.isclose(mass_residual, 0.0, abs_tol=1e-6)


def test_compressor_gamma_none_derives_from_fluid_and_matches_explicit_gamma():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1005.0 / (1005.0 - 287.05)  # what air.gamma(P, T) computes internally
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)
    N = 50.0 * math.sqrt(T_in) * 60.0
    mdot = 2.0 * (P_in / 1.0e5) / math.sqrt(T_in)
    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={
            "c1.in": (P_in, h_in),
            "c1.out": (3.0 * P_in, h_in),  # arbitrary guess, same for both
        },
    )

    explicit = Compressor(name="c1", map_path=_MAP_PATH, N=N, gamma=gamma)
    derived = Compressor(name="c1", map_path=_MAP_PATH, N=N, gamma=None)
    assert derived.gamma is None

    assert derived.residuals(state) == explicit.residuals(state)


def test_compressor_factor_overrides_defaults_to_file_values():
    comp_default = Compressor(name="c1", map_path=_MAP_PATH, N=3000.0, gamma=1.4)
    comp_explicit_none = Compressor(
        name="c1", map_path=_MAP_PATH, N=3000.0, gamma=1.4, factor_overrides=None
    )
    A, B = comp_default.map.mid_speed(), 1.0
    assert math.isclose(
        comp_default.map.pressure_ratio(A, B), comp_explicit_none.map.pressure_ratio(A, B)
    )


def test_compressor_factor_overrides_calibrates_pressure_ratio():
    comp_default = Compressor(name="c1", map_path=_MAP_PATH, N=3000.0, gamma=1.4)
    comp_calibrated = Compressor(
        name="c1",
        map_path=_MAP_PATH,
        N=3000.0,
        gamma=1.4,
        factor_overrides={"C_fact": 1.05},
    )
    A, B = comp_default.map.mid_speed(), 1.0
    assert math.isclose(
        comp_calibrated.map.pressure_ratio(A, B),
        1.05 * comp_default.map.pressure_ratio(A, B),
        rel_tol=1e-9,
    )


def test_compressor_report_metrics_reflects_map_lookup():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1.4
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)

    N = 70.0 * math.sqrt(T_in) * 60.0
    comp = Compressor(name="c1", map_path=_MAP_PATH, N=N, gamma=gamma)
    mdot = 2.0 * (P_in / 1.0e5) / math.sqrt(T_in)

    PR, eta_s = 4.0, 0.85  # map values at A=70, B=2 (see fixture)
    T_out_isentropic = T_in * PR ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_in * PR, T_out_isentropic)
    h_out = h_in + (h_out_isentropic - h_in) / eta_s
    P_out = PR * P_in

    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={
            "c1.in": (P_in, h_in),
            "c1.out": (P_out, h_out),
        },
    )
    metrics = comp.report_metrics(state)
    assert math.isclose(metrics["PR [-]"], PR, rel_tol=1e-9)
    assert math.isclose(metrics["eta_s [-]"], eta_s, rel_tol=1e-9)
    assert math.isclose(metrics["power [W]"], mdot * (h_out - h_in), rel_tol=1e-9)


def test_compressor_free_parameters_empty_when_n_given():
    comp = Compressor(name="c1", map_path=_MAP_PATH, N=3000.0, gamma=1.4)
    assert comp.free_parameters() == {}


def test_compressor_free_parameters_includes_n_when_n_omitted():
    comp = Compressor(name="c1", map_path=_MAP_PATH, gamma=1.4)
    params = comp.free_parameters()
    assert set(params.keys()) == {"N"}
    assert params["N"] > 0


def test_compressor_residuals_uses_n_from_state_param_when_n_omitted():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1.4
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)

    # A = 50 exactly, at T_in = 300 K -> N such that this A results.
    N = 50.0 * math.sqrt(T_in) * 60.0
    mdot = 2.0 * (P_in / 1.0e5) / math.sqrt(T_in)  # B = 2 exactly

    comp = Compressor(name="c1", map_path=_MAP_PATH, gamma=gamma)  # N left free
    PR, eta_s = 3.0, 0.8  # map values at A=50, B=2 (see fixture)
    T_out_isentropic = T_in * PR ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_in * PR, T_out_isentropic)
    h_out_guess = h_in + (h_out_isentropic - h_in) / eta_s

    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={
            "c1.in": (P_in, h_in),
            "c1.out": (PR * P_in, h_out_guess),
        },
        params={"c1.N": N},
    )
    momentum_residual, energy_residual, mass_residual = comp.residuals(state)
    assert math.isclose(momentum_residual, 0.0, abs_tol=1e-6)
    assert math.isclose(energy_residual, 0.0, abs_tol=1e-3)
    assert math.isclose(mass_residual, 0.0, abs_tol=1e-6)


def test_compressor_rejects_gamma_not_greater_than_one():
    try:
        Compressor(name="c1", map_path=_MAP_PATH, gamma=1.0)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_compressor_heat_path_defaults_to_none():
    comp = Compressor(name="c1", map_path=_MAP_PATH, N=3000.0, gamma=1.4)
    assert comp.heat_path is None


def test_compressor_heat_path_adds_q_over_mdot_to_energy_residual():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1.4
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)
    N = 50.0 * math.sqrt(T_in) * 60.0
    mdot = 2.0 * (P_in / 1.0e5) / math.sqrt(T_in)

    PR, eta_s = 3.0, 0.8  # map values at A=50, B=2 (see fixture)
    T_out_isentropic = T_in * PR ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_in * PR, T_out_isentropic)
    dh_actual = (h_out_isentropic - h_in) / eta_s
    h_out_guess = h_in + dh_actual

    state = _FakeState(
        fluid=air, mdot=mdot,
        node_values={"c1.in": (P_in, h_in), "c1.out": (PR * P_in, h_out_guess)},
    )

    baseline = Compressor(name="c1", map_path=_MAP_PATH, N=N, gamma=gamma)
    _, baseline_energy_residual, _ = baseline.residuals(state)

    Q = 3000.0
    lossy = Compressor(name="c1", map_path=_MAP_PATH, N=N, gamma=gamma)
    lossy.heat_path = _FakeHeatPath(Q_value=Q)
    _, lossy_energy_residual, _ = lossy.residuals(state)

    assert math.isclose(
        lossy_energy_residual - baseline_energy_residual, Q / mdot, rel_tol=1e-9
    )


def test_compressor_report_metrics_includes_q_loss():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)
    N = 50.0 * math.sqrt(T_in) * 60.0
    mdot = 2.0 * (P_in / 1.0e5) / math.sqrt(T_in)
    state = _FakeState(
        fluid=air, mdot=mdot, node_values={"c1.in": (P_in, h_in), "c1.out": (3.0 * P_in, h_in)},
    )

    comp = Compressor(name="c1", map_path=_MAP_PATH, N=N, gamma=1.4)
    assert comp.report_metrics(state)["Q_loss [W]"] == 0.0

    comp.heat_path = _FakeHeatPath(Q_value=777.0)
    assert comp.report_metrics(state)["Q_loss [W]"] == 777.0
