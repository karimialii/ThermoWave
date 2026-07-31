import math

import pytest

from thermowave.components.heat_transfer import Convection, ThermalMass
from thermowave.components.pipe import Pipe
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.core.settings import settings
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _FakeState:
    def __init__(self, fluid, mdot, node_values):
        self.fluid = fluid
        self._mdot = mdot
        self._node_values = node_values

    def node(self, name):
        return self._node_values[name]

    def fluid_at(self, name):
        return self.fluid

    def mdot(self, name):
        return self._mdot


def teardown_function(_):
    settings.friction_correlation = "haaland"


def test_f_given_directly_bypasses_roughness_calculation():
    pipe = Pipe(name="p1", L=5.0, D=0.2, f=0.02)
    assert pipe._friction_factor(rho=1.2, v=10.0) == 0.02


def test_missing_f_and_roughness_or_mu_raises():
    with pytest.raises(ValueError):
        Pipe(name="p1", L=5.0, D=0.2)
    with pytest.raises(ValueError):
        Pipe(name="p1", L=5.0, D=0.2, roughness=1e-5)
    with pytest.raises(ValueError):
        Pipe(name="p1", L=5.0, D=0.2, mu=1.8e-5)


def test_laminar_branch_matches_64_over_re():
    # Small D, small v, larger mu to keep Re << 2300.
    pipe = Pipe(name="p1", L=5.0, D=0.01, roughness=1e-6, mu=1.8e-3)
    rho, v = 1.2, 0.01
    Re = rho * v * pipe.D / pipe.mu
    assert Re < 2300.0
    f = pipe._friction_factor(rho, v)
    assert math.isclose(f, 64.0 / Re, rel_tol=1e-9)


def test_turbulent_haaland_branch_matches_hand_calc():
    pipe = Pipe(name="p1", L=5.0, D=0.1, roughness=4.5e-5, mu=1.8e-5)
    rho, v = 1.2, 20.0
    Re = rho * v * pipe.D / pipe.mu
    assert Re >= 2300.0
    rel_roughness = pipe.roughness / pipe.D
    x = -1.8 * math.log10((rel_roughness / 3.7) ** 1.11 + 6.9 / Re)
    expected = 1.0 / x**2
    assert math.isclose(pipe._friction_factor(rho, v), expected, rel_tol=1e-9)


def test_turbulent_colebrook_branch_close_to_but_not_identical_to_haaland():
    pipe = Pipe(name="p1", L=5.0, D=0.1, roughness=4.5e-5, mu=1.8e-5)
    rho, v = 1.2, 20.0

    settings.friction_correlation = "haaland"
    f_haaland = pipe._friction_factor(rho, v)

    settings.friction_correlation = "colebrook"
    f_colebrook = pipe._friction_factor(rho, v)

    assert not math.isclose(f_haaland, f_colebrook, rel_tol=1e-9)
    assert math.isclose(f_haaland, f_colebrook, rel_tol=0.02)


def test_settings_friction_correlation_switches_branch():
    pipe = Pipe(name="p1", L=5.0, D=0.1, roughness=4.5e-5, mu=1.8e-5)
    rho, v = 1.2, 20.0

    settings.friction_correlation = "haaland"
    f1 = pipe._friction_factor(rho, v)
    settings.friction_correlation = "colebrook"
    f2 = pipe._friction_factor(rho, v)
    assert f1 != f2


def test_invalid_friction_correlation_raises():
    pipe = Pipe(name="p1", L=5.0, D=0.1, roughness=4.5e-5, mu=1.8e-5)
    settings.friction_correlation = "not_a_real_correlation"
    with pytest.raises(ValueError):
        pipe._friction_factor(rho=1.2, v=20.0)


def test_residuals_use_computed_friction_factor():
    pipe = Pipe(name="p1", L=5.0, D=0.2, roughness=1e-5, mu=1.8e-5)
    area = math.pi * 0.2**2 / 4
    mdot = 1.0
    P_in, h_in = 200_000.0, 300_000.0
    state = _FakeState(
        fluid=AIR,
        mdot=mdot,
        node_values={"p1.in": (P_in, h_in), "p1.out": (190_000.0, h_in)},
    )
    rho = AIR.density_ph(P_in, h_in)
    v = mdot / (rho * area)
    f = pipe._friction_factor(rho, v)
    dp_expected = f * (pipe.L / pipe.D) * (rho * v**2 / 2)

    residuals = pipe.residuals(state)
    momentum_residual = residuals[0]
    assert math.isclose(
        momentum_residual, P_in - 190_000.0 - dp_expected, rel_tol=1e-9
    )


def test_heat_path_none_matches_pre_existing_behavior():
    # heat_path defaults to None; heat_loss_watts(None, ...) short-circuits
    # to 0.0 without ever touching state, so this must reproduce the exact
    # heat_loss-only residual every other test in this file already checks.
    pipe = Pipe(name="p1", L=5.0, D=0.2, f=0.02, heat_loss=-1000.0)
    state = _FakeState(
        fluid=AIR, mdot=1.0,
        node_values={"p1.in": (200_000.0, 300_000.0), "p1.out": (190_000.0, 300_000.0)},
    )
    energy_residual = pipe.residuals(state)[1]
    assert math.isclose(energy_residual, 0.0 - (-1000.0 / 1.0), rel_tol=1e-9)


def test_pipe_accepts_heat_path_via_add_heat_path():
    # A Pipe attached to a live Convection/ThermalMass heat path (rather
    # than a fixed heat_loss scalar) -- the mechanism the combustion-chamber
    # liner guide relies on to couple each discretized station to its own
    # wall ring.
    src = Source(name="src", P=200_000.0, T=400.0, mdot=1.0)
    pipe = Pipe(name="p1", L=1.0, D=0.1, f=0.02, n_elem=1)
    snk = Sink(name="snk")
    wall = ThermalMass(name="wall", thermal_capacitance=500.0, T0=300.0)

    network = Network(fluid=AIR)
    for c in (src, pipe, snk):
        network.add_component(c)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")
    network.add_heat_path(Convection(name="conv", a=(pipe, "in"), b=wall, h=50.0, A=0.3))
    # A second path (wall -> fixed ambient) so equilibrium isn't the trivial
    # zero-net-Q case of a mass with only one heat source -- otherwise the
    # wall simply settles at the gas temperature and no heat actually flows.
    network.add_heat_path(Convection(name="conv_ambient", a=wall, b=250.0, h=5.0, A=0.3))

    assert network.check_wiring() == []
    result = network.solve(tol=1e-9, max_iter=200, verbose=False, progress=False)
    assert result.converged

    state = result.state()
    T_wall = state.param("wall.T")
    T_pipe_in = state.fluid_at("p1.in").temperature_ph(*state.node("p1.in"))
    # Steady state: the wall's net heat is zero across BOTH its paths, so it
    # settles strictly between the hot gas it convects with and the fixed
    # ambient it also loses heat to -- not pinned to either endpoint.
    assert 250.0 < T_wall < T_pipe_in

    # And the pipe's own energy balance actually lost heat to get there:
    # h_out should be below h_in by the convected Q/mdot.
    from thermowave.components.heat_transfer import heat_loss_watts
    Q_loss = heat_loss_watts(pipe.heat_path, state)
    assert Q_loss > 0.0
    P_out, h_out = state.node("p1.out")
    P_in, h_in = state.node("p1.in")
    assert math.isclose(h_in - h_out, Q_loss / 1.0, rel_tol=1e-6)
