import math

from thermowave.components.pipe import Pipe
from thermowave.core.constants import STANDARD_GRAVITY_M_S2
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


def test_default_elevation_is_zero_no_gravity_term():
    # No z_in/z_out given -- momentum residual must match the pre-existing
    # friction-only behavior exactly (regression: every caller that never
    # heard of elevation shouldn't see any change).
    pipe = Pipe(name="p1", L=5.0, D=0.2, f=0.02)
    assert pipe.z_in == 0.0
    assert pipe.z_out == 0.0

    area = math.pi * 0.2**2 / 4
    mdot = 1.0
    P_in, h_in = 200_000.0, 300_000.0
    state = _FakeState(
        fluid=AIR, mdot=mdot,
        node_values={"p1.in": (P_in, h_in), "p1.out": (190_000.0, h_in)},
    )
    rho = AIR.density_ph(P_in, h_in)
    v = mdot / (rho * area)
    dp_friction = 0.02 * (pipe.L / pipe.D) * (rho * v**2 / 2)

    momentum_residual = pipe.residuals(state)[0]
    assert math.isclose(momentum_residual, P_in - 190_000.0 - dp_friction, rel_tol=1e-9)


def test_gravity_term_added_for_rising_elevation():
    pipe = Pipe(name="p1", L=5.0, D=0.2, f=0.02, z_in=0.0, z_out=10.0)
    area = math.pi * 0.2**2 / 4
    mdot = 1.0
    P_in, h_in = 200_000.0, 300_000.0
    P_out = 190_000.0
    state = _FakeState(
        fluid=AIR, mdot=mdot,
        node_values={"p1.in": (P_in, h_in), "p1.out": (P_out, h_in)},
    )
    rho = AIR.density_ph(P_in, h_in)
    v = mdot / (rho * area)
    dp_friction = 0.02 * (pipe.L / pipe.D) * (rho * v**2 / 2)
    dp_gravity = rho * STANDARD_GRAVITY_M_S2 * (10.0 - 0.0)

    momentum_residual = pipe.residuals(state)[0]
    expected = P_in - P_out - dp_friction - dp_gravity
    assert math.isclose(momentum_residual, expected, rel_tol=1e-9)
    # Rising elevation must make the residual smaller than the no-gravity
    # case (an extra positive term is being subtracted).
    assert momentum_residual < P_in - P_out - dp_friction


def test_gravity_term_flips_sign_for_falling_elevation():
    pipe = Pipe(name="p1", L=5.0, D=0.2, f=0.02, z_in=10.0, z_out=0.0)
    area = math.pi * 0.2**2 / 4
    mdot = 1.0
    P_in, h_in = 200_000.0, 300_000.0
    P_out = 190_000.0
    state = _FakeState(
        fluid=AIR, mdot=mdot,
        node_values={"p1.in": (P_in, h_in), "p1.out": (P_out, h_in)},
    )
    rho = AIR.density_ph(P_in, h_in)
    v = mdot / (rho * area)
    dp_friction = 0.02 * (pipe.L / pipe.D) * (rho * v**2 / 2)
    dp_gravity = rho * STANDARD_GRAVITY_M_S2 * (0.0 - 10.0)  # negative

    momentum_residual = pipe.residuals(state)[0]
    expected = P_in - P_out - dp_friction - dp_gravity
    assert math.isclose(momentum_residual, expected, rel_tol=1e-9)
    assert momentum_residual > P_in - P_out - dp_friction


def test_gravity_term_split_evenly_across_elements():
    # n_elem=2 over a 20 m rise: each element should see a 10 m rise, not
    # the whole 20 m twice -- same "split evenly" contract elem_L already
    # follows for length.
    pipe = Pipe(name="p1", L=4.0, D=0.2, f=0.02, n_elem=2, z_in=0.0, z_out=20.0)
    area = math.pi * 0.2**2 / 4
    mdot = 1.0
    P0, h = 200_000.0, 300_000.0
    state = _FakeState(
        fluid=AIR, mdot=mdot,
        node_values={
            "p1.in": (P0, h),
            "p1__mid1": (P0 - 100.0, h),
            "p1.out": (P0 - 200.0, h),
        },
    )
    rho = AIR.density_ph(P0, h)
    v = mdot / (rho * area)
    dp_friction_elem = 0.02 * (2.0 / 0.2) * (rho * v**2 / 2)
    dp_gravity_elem = rho * STANDARD_GRAVITY_M_S2 * 10.0

    residuals = pipe.residuals(state)
    elem1_momentum = residuals[0]
    expected = P0 - (P0 - 100.0) - dp_friction_elem - dp_gravity_elem
    assert math.isclose(elem1_momentum, expected, rel_tol=1e-9)
