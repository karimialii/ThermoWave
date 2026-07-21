import math
import math as _math

from thermowave.components.pipe import Pipe
from thermowave.components.simple_compressor import SimpleCompressor
from thermowave.components.simple_turbine import SimpleTurbine
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.valve import Valve
from thermowave.core.settings import settings
from thermowave.fluids.ideal_gas import IdealGasFluid
from thermowave.fluids.ideal_gas import IdealGasFluid as _IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def test_source_ports_returns_outlet_port_derived_from_name():
    src = Source(name="src1", P=101325.0, T=300.0, mdot=1.0)
    assert src.ports() == {"out": "src1.out"}


def test_source_residuals_is_empty():
    src = Source(name="src1", P=101325.0, T=300.0, mdot=1.0)
    assert src.residuals(state=None) == []


def test_source_fixed_node_values_uses_si_units_by_default():
    src = Source(name="src1", P=101325.0, T=300.0, mdot=1.0)
    fixed = src.fixed_node_values(AIR)
    P, h = fixed["src1.out"]
    assert math.isclose(P, 101325.0, rel_tol=1e-9)
    assert math.isclose(h, 1005.0 * 300.0, rel_tol=1e-9)


def test_source_fixed_node_values_converts_non_si_units(reset_settings=None):
    settings.pressure_unit = "bar"
    settings.temperature_unit = "C"
    try:
        src = Source(name="src1", P=1.013, T=26.85, mdot=1.0)
        fixed = src.fixed_node_values(AIR)
        P, h = fixed["src1.out"]
        assert math.isclose(P, 101300.0, rel_tol=1e-9)
        assert math.isclose(h, 1005.0 * 300.0, rel_tol=1e-6)
    finally:
        settings.pressure_unit = "Pa"
        settings.temperature_unit = "K"


def test_source_fixed_node_mdot_returns_configured_value():
    src = Source(name="src1", P=101325.0, T=300.0, mdot=2.5)
    assert src.fixed_node_mdot() == {"src1.out": 2.5}


def test_sink_ports_returns_inlet_port_derived_from_name():
    snk = Sink(name="snk1")
    assert snk.ports() == {"in": "snk1.in"}


def test_sink_residuals_is_empty():
    snk = Sink(name="snk1")
    assert snk.residuals(state=None) == []


def test_sink_has_no_fixed_node_values_or_mdot():
    snk = Sink(name="snk1")
    assert snk.fixed_node_values(AIR) == {}
    assert snk.fixed_node_mdot() == {}


class _FakeState:
    """Minimal stand-in for core.network.NetworkState.

    mdot may be a single float (broadcast to every node) or a dict keyed by
    node name, for tests that need distinct flow rates on different ports.
    """

    def __init__(self, fluid, mdot, node_values: dict[str, tuple[float, float]]):
        self.fluid = fluid
        self._mdot = mdot
        self._node_values = node_values

    def node(self, name: str) -> tuple[float, float]:
        return self._node_values[name]

    def fluid_at(self, name: str):
        return self.fluid

    def mdot(self, name: str) -> float:
        if isinstance(self._mdot, dict):
            return self._mdot[name]
        return self._mdot


def test_pipe_ports_returns_inlet_and_outlet_derived_from_name():
    pipe = Pipe(name="p1", L=5.0, D=0.2, f=0.02, n_elem=1)
    assert pipe.ports() == {"in": "p1.in", "out": "p1.out"}


def test_pipe_internal_nodes_empty_for_single_element():
    pipe = Pipe(name="p1", L=5.0, D=0.2, f=0.02, n_elem=1)
    assert pipe.internal_nodes() == []


def test_pipe_internal_nodes_for_multiple_elements():
    pipe = Pipe(name="p1", L=5.0, D=0.2, f=0.02, n_elem=3)
    assert pipe.internal_nodes() == ["p1__mid1", "p1__mid2"]


def test_pipe_residuals_single_element_matches_darcy_weisbach_hand_calc():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    L, D, f = 5.0, 0.2, 0.02
    mdot = 1.0
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)

    rho = air.density_ph(P_in, h_in)
    area = _math.pi * D**2 / 4
    v = mdot / (rho * area)
    expected_dp = f * (L / D) * (rho * v**2 / 2)
    P_out_guess = P_in - expected_dp
    h_out_guess = h_in  # adiabatic

    pipe = Pipe(name="p1", L=L, D=D, f=f, n_elem=1)
    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={
            "p1.in": (P_in, h_in),
            "p1.out": (P_out_guess, h_out_guess),
        },
    )
    residuals = pipe.residuals(state)
    assert len(residuals) == 3
    momentum_residual, energy_residual, mass_residual = residuals
    assert _math.isclose(momentum_residual, 0.0, abs_tol=1e-6)
    assert _math.isclose(energy_residual, 0.0, abs_tol=1e-6)
    assert _math.isclose(mass_residual, 0.0, abs_tol=1e-6)


def test_pipe_residuals_with_heat_loss_reflects_energy_balance():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    L, D, f = 5.0, 0.2, 0.02
    mdot = 1.0
    Q = 1000.0  # W, total heat loss
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)
    h_out_guess = h_in - Q / mdot  # exact energy balance
    P_out_guess = P_in  # momentum residual will be nonzero; only checking energy term

    pipe = Pipe(name="p1", L=L, D=D, f=f, n_elem=1, heat_loss=Q)
    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={
            "p1.in": (P_in, h_in),
            "p1.out": (P_out_guess, h_out_guess),
        },
    )
    _momentum_residual, energy_residual, _mass_residual = pipe.residuals(state)
    assert _math.isclose(energy_residual, 0.0, abs_tol=1e-6)


def test_pipe_residuals_multi_element_returns_two_residuals_per_element():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    pipe = Pipe(name="p1", L=6.0, D=0.2, f=0.02, n_elem=3)
    P_in, h_in = 101325.0, air.enthalpy_pt(101325.0, 300.0)
    state = _FakeState(
        fluid=air,
        mdot=1.0,
        node_values={
            "p1.in": (P_in, h_in),
            "p1__mid1": (P_in - 100.0, h_in),
            "p1__mid2": (P_in - 200.0, h_in),
            "p1.out": (P_in - 300.0, h_in),
        },
    )
    residuals = pipe.residuals(state)
    assert len(residuals) == 7  # 2 residuals (momentum, energy) * 3 elements + 1 mass


def test_valve_ports_returns_inlet_and_outlet_derived_from_name():
    valve = Valve(name="v1", D=0.1, K=5.0)
    assert valve.ports() == {"in": "v1.in", "out": "v1.out"}


def test_valve_residuals_matches_k_factor_hand_calc_fully_open():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    D, K = 0.1, 5.0
    mdot = 0.5
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)

    area = _math.pi * D**2 / 4
    rho = air.density_ph(P_in, h_in)
    v = mdot / (rho * area)
    expected_dp = K * (rho * v**2 / 2)
    P_out_guess = P_in - expected_dp
    h_out_guess = h_in  # isenthalpic

    valve = Valve(name="v1", D=D, K=K)
    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={
            "v1.in": (P_in, h_in),
            "v1.out": (P_out_guess, h_out_guess),
        },
    )
    momentum_residual, energy_residual, mass_residual = valve.residuals(state)
    assert _math.isclose(momentum_residual, 0.0, abs_tol=1e-6)
    assert _math.isclose(energy_residual, 0.0, abs_tol=1e-6)
    assert _math.isclose(mass_residual, 0.0, abs_tol=1e-6)


def test_valve_residuals_scales_k_by_inverse_square_of_opening():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    D, K, opening = 0.1, 5.0, 0.5
    mdot = 0.5
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)

    area = _math.pi * D**2 / 4
    rho = air.density_ph(P_in, h_in)
    v = mdot / (rho * area)
    k_eff = K / opening**2
    expected_dp = k_eff * (rho * v**2 / 2)
    P_out_guess = P_in - expected_dp

    valve = Valve(name="v1", D=D, K=K, opening=opening)
    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={
            "v1.in": (P_in, h_in),
            "v1.out": (P_out_guess, h_in),
        },
    )
    momentum_residual, _energy_residual, _mass_residual = valve.residuals(state)
    assert _math.isclose(momentum_residual, 0.0, abs_tol=1e-6)


def test_valve_residuals_energy_residual_reflects_enthalpy_mismatch():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)

    valve = Valve(name="v1", D=0.1, K=5.0)
    state = _FakeState(
        fluid=air,
        mdot=0.5,
        node_values={
            "v1.in": (P_in, h_in),
            "v1.out": (P_in, h_in - 500.0),  # wrong: valve should be isenthalpic
        },
    )
    _momentum_residual, energy_residual, _mass_residual = valve.residuals(state)
    assert _math.isclose(energy_residual, 500.0, abs_tol=1e-6)


def test_simple_compressor_ports_returns_inlet_and_outlet_derived_from_name():
    comp = SimpleCompressor(name="c1", PR=3.0, eta_s=0.8, gamma=1.4)
    assert comp.ports() == {"in": "c1.in", "out": "c1.out"}


def test_simple_compressor_residuals_matches_isentropic_efficiency_hand_calc():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1005.0 / (1005.0 - 287.05)  # matches air's cp/R exactly
    PR, eta_s = 3.0, 0.8
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)

    T_out_isentropic = T_in * PR ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_in * PR, T_out_isentropic)
    dh_actual = (h_out_isentropic - h_in) / eta_s
    P_out_guess = PR * P_in
    h_out_guess = h_in + dh_actual

    comp = SimpleCompressor(name="c1", PR=PR, eta_s=eta_s, gamma=gamma)
    state = _FakeState(
        fluid=air,
        mdot=1.0,
        node_values={
            "c1.in": (P_in, h_in),
            "c1.out": (P_out_guess, h_out_guess),
        },
    )
    momentum_residual, energy_residual, mass_residual = comp.residuals(state)
    assert _math.isclose(momentum_residual, 0.0, abs_tol=1e-6)
    assert _math.isclose(energy_residual, 0.0, abs_tol=1e-3)
    assert _math.isclose(mass_residual, 0.0, abs_tol=1e-6)


def test_simple_compressor_gamma_none_derives_from_fluid_and_matches_explicit_gamma():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1005.0 / (1005.0 - 287.05)  # what air.gamma(P, T) computes internally
    PR, eta_s = 3.0, 0.8
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)
    state = _FakeState(
        fluid=air,
        mdot=1.0,
        node_values={
            "c1.in": (P_in, h_in),
            "c1.out": (PR * P_in, h_in),  # arbitrary guess, same for both
        },
    )

    explicit = SimpleCompressor(name="c1", PR=PR, eta_s=eta_s, gamma=gamma)
    derived = SimpleCompressor(name="c1", PR=PR, eta_s=eta_s, gamma=None)
    assert derived.gamma is None  # confirms it wasn't silently defaulted at construction

    assert derived.residuals(state) == explicit.residuals(state)


def test_simple_compressor_actual_enthalpy_rise_exceeds_isentropic_when_inefficient():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1.4
    PR = 2.0
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)

    T_out_isentropic = T_in * PR ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_in * PR, T_out_isentropic)
    dh_isentropic = h_out_isentropic - h_in

    comp = SimpleCompressor(name="c1", PR=PR, eta_s=0.75, gamma=gamma)
    state = _FakeState(
        fluid=air,
        mdot=1.0,
        # Guess the isentropic (not actual) outlet enthalpy: with eta_s < 1
        # the true actual rise must be larger, so this guess is too low and
        # the energy residual (guess - required) should be negative.
        node_values={
            "c1.in": (P_in, h_in),
            "c1.out": (PR * P_in, h_out_isentropic),
        },
    )
    _momentum_residual, energy_residual, _mass_residual = comp.residuals(state)
    assert energy_residual < 0
    assert _math.isclose(energy_residual, dh_isentropic * (1 - 1 / 0.75), abs_tol=1e-3)


class _FakeHeatPath:
    def __init__(self, Q_value):
        self._Q_value = Q_value

    def Q(self, state):
        return self._Q_value


def test_simple_compressor_heat_path_defaults_to_none():
    comp = SimpleCompressor(name="c1", PR=3.0, eta_s=0.8, gamma=1.4)
    assert comp.heat_path is None


def test_simple_compressor_heat_path_adds_q_over_mdot_to_energy_residual():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1.4
    PR, eta_s = 3.0, 0.8
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)
    T_out_isentropic = T_in * PR ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_in * PR, T_out_isentropic)
    h_out_guess = h_in + (h_out_isentropic - h_in) / eta_s
    mdot = 1.0

    state = _FakeState(
        fluid=air, mdot=mdot,
        node_values={"c1.in": (P_in, h_in), "c1.out": (PR * P_in, h_out_guess)},
    )

    baseline = SimpleCompressor(name="c1", PR=PR, eta_s=eta_s, gamma=gamma)
    _, baseline_energy_residual, _ = baseline.residuals(state)

    Q = 2000.0
    lossy = SimpleCompressor(name="c1", PR=PR, eta_s=eta_s, gamma=gamma)
    lossy.heat_path = _FakeHeatPath(Q_value=Q)
    _, lossy_energy_residual, _ = lossy.residuals(state)

    assert _math.isclose(
        lossy_energy_residual - baseline_energy_residual, Q / mdot, rel_tol=1e-9
    )


def test_simple_turbine_ports_returns_inlet_and_outlet_derived_from_name():
    turb = SimpleTurbine(name="t1", PR=2.0, eta_s=0.85, gamma=1.4)
    assert turb.ports() == {"in": "t1.in", "out": "t1.out"}


def test_simple_turbine_residuals_matches_isentropic_efficiency_hand_calc():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1005.0 / (1005.0 - 287.05)  # matches air's cp/R exactly
    PR, eta_s = 2.0, 0.85
    P_in, T_in = 400000.0, 1000.0
    h_in = air.enthalpy_pt(P_in, T_in)
    P_out = P_in / PR

    T_out_isentropic = T_in * (1.0 / PR) ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_out, T_out_isentropic)
    dh_actual = eta_s * (h_in - h_out_isentropic)
    h_out_guess = h_in - dh_actual

    turb = SimpleTurbine(name="t1", PR=PR, eta_s=eta_s, gamma=gamma)
    state = _FakeState(
        fluid=air,
        mdot=1.0,
        node_values={
            "t1.in": (P_in, h_in),
            "t1.out": (P_out, h_out_guess),
        },
    )
    momentum_residual, energy_residual, mass_residual = turb.residuals(state)
    assert _math.isclose(momentum_residual, 0.0, abs_tol=1e-6)
    assert _math.isclose(energy_residual, 0.0, abs_tol=1e-3)
    assert _math.isclose(mass_residual, 0.0, abs_tol=1e-6)


def test_simple_turbine_gamma_none_derives_from_fluid_and_matches_explicit_gamma():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1005.0 / (1005.0 - 287.05)
    PR, eta_s = 2.0, 0.85
    P_in, T_in = 400000.0, 1000.0
    h_in = air.enthalpy_pt(P_in, T_in)
    P_out = P_in / PR
    state = _FakeState(
        fluid=air,
        mdot=1.0,
        node_values={
            "t1.in": (P_in, h_in),
            "t1.out": (P_out, h_in),  # arbitrary guess, same for both
        },
    )

    explicit = SimpleTurbine(name="t1", PR=PR, eta_s=eta_s, gamma=gamma)
    derived = SimpleTurbine(name="t1", PR=PR, eta_s=eta_s, gamma=None)
    assert derived.gamma is None

    assert derived.residuals(state) == explicit.residuals(state)


def test_simple_turbine_actual_enthalpy_drop_is_less_than_isentropic_when_inefficient():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1.4
    PR = 2.0
    P_in, T_in = 400000.0, 1000.0
    h_in = air.enthalpy_pt(P_in, T_in)
    P_out = P_in / PR

    T_out_isentropic = T_in * (1.0 / PR) ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_out, T_out_isentropic)
    dh_isentropic = h_in - h_out_isentropic

    turb = SimpleTurbine(name="t1", PR=PR, eta_s=0.75, gamma=gamma)
    state = _FakeState(
        fluid=air,
        mdot=1.0,
        # Guess the isentropic (not actual) outlet enthalpy: with eta_s < 1
        # the true actual drop must be smaller, so this guess is too low and
        # the energy residual (guess - required) should be negative.
        node_values={
            "t1.in": (P_in, h_in),
            "t1.out": (P_out, h_out_isentropic),
        },
    )
    _momentum_residual, energy_residual, _mass_residual = turb.residuals(state)
    assert energy_residual < 0
    assert _math.isclose(energy_residual, -dh_isentropic * (1 - 0.75), abs_tol=1e-3)


def test_simple_turbine_heat_path_defaults_to_none():
    turb = SimpleTurbine(name="t1", PR=2.0, eta_s=0.85, gamma=1.4)
    assert turb.heat_path is None


def test_simple_turbine_heat_path_adds_q_over_mdot_to_energy_residual():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    gamma = 1.4
    PR, eta_s = 2.0, 0.85
    P_in, T_in = 400000.0, 1000.0
    h_in = air.enthalpy_pt(P_in, T_in)
    P_out = P_in / PR
    T_out_isentropic = T_in * (1.0 / PR) ** ((gamma - 1.0) / gamma)
    h_out_isentropic = air.enthalpy_pt(P_out, T_out_isentropic)
    h_out_guess = h_in - eta_s * (h_in - h_out_isentropic)
    mdot = 1.0

    state = _FakeState(
        fluid=air, mdot=mdot, node_values={"t1.in": (P_in, h_in), "t1.out": (P_out, h_out_guess)},
    )

    baseline = SimpleTurbine(name="t1", PR=PR, eta_s=eta_s, gamma=gamma)
    _, baseline_energy_residual, _ = baseline.residuals(state)

    Q = 1500.0
    lossy = SimpleTurbine(name="t1", PR=PR, eta_s=eta_s, gamma=gamma)
    lossy.heat_path = _FakeHeatPath(Q_value=Q)
    _, lossy_energy_residual, _ = lossy.residuals(state)

    assert _math.isclose(
        lossy_energy_residual - baseline_energy_residual, Q / mdot, rel_tol=1e-9
    )


def test_simple_compressor_rejects_invalid_eta_s():
    try:
        SimpleCompressor(name="c1", PR=3.0, eta_s=1.5, gamma=1.4)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_simple_compressor_rejects_non_positive_pr():
    try:
        SimpleCompressor(name="c1", PR=0.0, eta_s=0.8, gamma=1.4)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_simple_compressor_rejects_gamma_not_greater_than_one():
    try:
        SimpleCompressor(name="c1", PR=3.0, eta_s=0.8, gamma=1.0)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_simple_turbine_rejects_invalid_eta_s():
    try:
        SimpleTurbine(name="t1", PR=2.0, eta_s=0.0, gamma=1.4)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_simple_turbine_rejects_non_positive_pr():
    try:
        SimpleTurbine(name="t1", PR=-1.0, eta_s=0.85, gamma=1.4)
        assert False, "expected ValueError"
    except ValueError:
        pass
