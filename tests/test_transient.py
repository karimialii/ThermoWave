import copy
import math

import pytest

from thermowave.components.compressor import Compressor
from thermowave.components.pipe import Pipe
from thermowave.components.shaft import Shaft
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)
GAMMA = 1005.0 / (1005.0 - 287.05)


def _build_dynamic_turboshaft(N0: float = 60000.0, inertia: float = 0.05):
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    comp = Compressor(name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None)
    heater = Pipe(name="heater", L=1.0, D=0.1, f=0.0, n_elem=1, heat_loss=-300000.0)
    turb = Turbine(name="turb", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=None)
    shaft = Shaft(
        name="shaft", components=[comp, turb], signs=[-1.0, 1.0],
        efficiency=0.98, inertia=inertia, dynamic=True, N0=N0,
    )
    snk = Sink(name="snk")

    network = Network(fluid=AIR)
    for component in (src, comp, heater, turb, shaft, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", heater, "in")
    network.connect(heater, "out", turb, "in")
    network.connect(turb, "out", snk, "in")
    return network, shaft


def test_solve_transient_rejects_network_with_nothing_time_varying():
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    pipe = Pipe(name="pipe", L=1.0, D=0.1, f=0.02)
    snk = Sink(name="snk")
    network = Network(fluid=AIR)
    for component in (src, pipe, snk):
        network.add_component(component)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")

    with pytest.raises(ValueError, match="nothing to evolve"):
        network.solve_transient(duration=1.0, dt=0.1)


def test_solve_transient_default_initial_condition_is_already_at_equilibrium():
    # No initial= given -> solve_transient() runs an ordinary steady-state
    # solve first, which (since Shaft's speed closes via derivative == 0)
    # is already the torque-balance equilibrium: speed shouldn't move at all
    # across subsequent transient steps.
    network, shaft = _build_dynamic_turboshaft(N0=60000.0)
    history = network.solve_transient(
        duration=0.5, dt=0.1, tol=1e-8, max_iter=400, damping=0.3,
    )
    speeds = history.diff_history["shaft.N"]
    assert len(speeds) == len(history.times) == len(history.steps)
    for N in speeds[1:]:
        assert math.isclose(N, speeds[0], rel_tol=1e-4)


def test_solve_transient_history_uses_initial_as_t0():
    network, shaft = _build_dynamic_turboshaft(N0=60000.0)
    initial = network.solve(tol=1e-8, max_iter=400, damping=0.3)

    history = network.solve_transient(
        duration=0.2, dt=0.05, initial=initial, tol=1e-8, max_iter=400, damping=0.3,
    )
    assert history.steps[0] is initial
    assert history.times[0] == 0.0
    assert math.isclose(
        history.diff_history["shaft.N"][0], initial.params["shaft.N"]
    )
    assert math.isclose(history.times[-1], 0.2, rel_tol=1e-9)


def test_solve_transient_off_equilibrium_initial_condition_satisfies_backward_euler():
    network, shaft = _build_dynamic_turboshaft(N0=60000.0)
    equilibrium = network.solve(tol=1e-8, max_iter=400, damping=0.3)

    # Craft a deliberately off-equilibrium t=0 state by starting the rotor
    # well below the torque-balance speed found above.
    off_equilibrium = copy.copy(equilibrium)
    off_equilibrium.params = dict(equilibrium.params)
    off_equilibrium.params["shaft.N"] = 50000.0

    dt = 0.05
    history = network.solve_transient(
        duration=dt, dt=dt, initial=off_equilibrium, tol=1e-8, max_iter=400, damping=0.3,
    )
    N0, N1 = history.diff_history["shaft.N"]
    assert math.isclose(N0, 50000.0)
    assert N1 > N0  # below equilibrium -> net positive torque -> speeds up

    state = NetworkState(
        fluid=history.steps[1].fluid, node_P=history.steps[1].node_P,
        node_h=history.steps[1].node_h, node_mdot=history.steps[1].node_mdot,
        params=history.steps[1].params,
    )
    rate = shaft.state_derivative(state)["N"]
    assert math.isclose((N1 - N0) / dt, rate, rel_tol=1e-4)


def test_solve_transient_diff_history_matches_step_params():
    network, shaft = _build_dynamic_turboshaft(N0=60000.0)
    history = network.solve_transient(
        duration=0.15, dt=0.05, tol=1e-8, max_iter=400, damping=0.3,
    )
    for N, step in zip(history.diff_history["shaft.N"], history.steps):
        assert math.isclose(N, step.params["shaft.N"])
