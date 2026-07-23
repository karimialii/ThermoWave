import copy
import io
import math

import pytest

from thermowave.components.compressor import Compressor
from thermowave.components.pid_controller import PIDController
from thermowave.components.pipe import Pipe
from thermowave.components.sensor import Sensor
from thermowave.components.shaft import Shaft
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.exceptions import ConvergenceError
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


def _off_equilibrium_initial(network, shaft, N0_target: float):
    equilibrium = network.solve(tol=1e-8, max_iter=400, damping=0.3)
    off_equilibrium = copy.copy(equilibrium)
    off_equilibrium.params = dict(equilibrium.params)
    off_equilibrium.params["shaft.N"] = N0_target
    return off_equilibrium


def test_solve_transient_adaptive_matches_fine_fixed_step_reference():
    # A loose-tolerance adaptive run should land close to a fine fixed-step
    # reference over the same off-equilibrium spin-up transient.
    network, shaft = _build_dynamic_turboshaft(N0=60000.0)
    initial = _off_equilibrium_initial(network, shaft, N0_target=50000.0)

    adaptive_history = network.solve_transient(
        duration=0.5, dt=0.01, initial=initial, tol=1e-8, max_iter=400, damping=0.3,
        adaptive=True, rtol=1e-4, atol=1.0,
    )

    network2, shaft2 = _build_dynamic_turboshaft(N0=60000.0)
    initial2 = _off_equilibrium_initial(network2, shaft2, N0_target=50000.0)
    fine_history = network2.solve_transient(
        duration=0.5, dt=0.002, initial=initial2, tol=1e-8, max_iter=400, damping=0.3,
    )

    assert math.isclose(
        adaptive_history.diff_history["shaft.N"][-1],
        fine_history.diff_history["shaft.N"][-1],
        rel_tol=1e-3,
    )
    assert math.isclose(adaptive_history.times[-1], 0.5, rel_tol=1e-9)
    # The point of adaptive stepping: fewer steps than the fine reference for
    # this smoothly-relaxing transient, since it can grow dt as the
    # transient settles.
    assert len(adaptive_history.times) < len(fine_history.times)


def test_solve_transient_adaptive_requires_a_differential_state():
    # A PID-only network (no dynamic Shaft/Tank) has nothing for
    # step-doubling's error estimate to compare against.
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    comp = Compressor(name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None)
    sensor = Sensor(name="s1")
    snk = Sink(name="snk")
    pid = PIDController(
        name="pid", sensor=sensor, quantity="T [K]", component=comp, free_param="N",
        setpoint=420.0, Kp=60.0, Ki=50.0, Kd=0.0, output0=60000.0,
    )
    network = Network(fluid=AIR)
    for component in (src, comp, sensor, pid, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", sensor, "tap")
    network.connect(comp, "out", snk, "in")

    with pytest.raises(ValueError, match="differential state"):
        network.solve_transient(duration=1.0, dt=0.1, adaptive=True)


def test_solve_transient_adaptive_respects_dt_max():
    network, shaft = _build_dynamic_turboshaft(N0=60000.0)
    initial = _off_equilibrium_initial(network, shaft, N0_target=50000.0)

    dt_max = 0.02
    history = network.solve_transient(
        duration=0.3, dt=0.005, initial=initial, tol=1e-8, max_iter=400, damping=0.3,
        adaptive=True, rtol=1e-4, atol=1.0, dt_max=dt_max,
    )
    step_sizes = [b - a for a, b in zip(history.times[:-1], history.times[1:])]
    assert all(h <= dt_max + 1e-9 for h in step_sizes)


def test_solve_transient_adaptive_gives_up_after_max_step_shrinks_if_unmeetable():
    # An unreasonably tight rtol/atol with dt_min pinned well above what's
    # needed to meet it should raise ConvergenceError rather than looping
    # forever trying to shrink past dt_min.
    network, shaft = _build_dynamic_turboshaft(N0=60000.0)
    initial = _off_equilibrium_initial(network, shaft, N0_target=30000.0)

    with pytest.raises(ConvergenceError, match="rejected"):
        network.solve_transient(
            duration=0.5, dt=0.05, initial=initial, tol=1e-8, max_iter=400, damping=0.3,
            adaptive=True, rtol=1e-12, atol=1e-12, dt_min=0.04, max_step_shrinks=3,
        )


def test_solve_transient_adaptive_calls_step_once_per_accepted_step_with_accepted_dt():
    # step() (a PIDController's finite-response update) must never fire for
    # a rejected/trial step, and must use that step's real (post-control)
    # dt, not the pre-adaptive initial guess. The PID loop here is wired as
    # a second, independent branch alongside the dynamic turboshaft (its own
    # Source/Compressor/Sink) so its free parameter doesn't collide with the
    # shaft's own speed-tie residuals.
    network, shaft = _build_dynamic_turboshaft(N0=50000.0)

    src2 = Source(name="src2", P=101325.0, T=288.15, mdot=0.5)
    comp2 = Compressor(name="comp2", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None)
    sensor = Sensor(name="s1")
    snk2 = Sink(name="snk2")
    pid = PIDController(
        name="pid", sensor=sensor, quantity="T [K]", component=comp2, free_param="N",
        setpoint=420.0, Kp=60.0, Ki=50.0, Kd=0.0, output0=60000.0,
    )
    for component in (src2, comp2, sensor, pid, snk2):
        network.add_component(component)
    network.connect(src2, "out", comp2, "in")
    network.connect(comp2, "out", sensor, "tap")
    network.connect(comp2, "out", snk2, "in")

    initial = _off_equilibrium_initial(network, shaft, N0_target=45000.0)

    recorded_dts = []
    original_step = PIDController.step

    def _spy_step(self, state, dt):
        recorded_dts.append(dt)
        return original_step(self, state, dt)

    PIDController.step = _spy_step
    try:
        history = network.solve_transient(
            duration=0.3, dt=0.02, initial=initial, tol=1e-8, max_iter=400, damping=0.3,
            adaptive=True, rtol=1e-3, atol=1.0,
        )
    finally:
        PIDController.step = original_step

    actual_dts = [b - a for a, b in zip(history.times[:-1], history.times[1:])]
    assert len(recorded_dts) == len(actual_dts)
    for recorded, actual in zip(recorded_dts, actual_dts):
        assert math.isclose(recorded, actual, rel_tol=1e-9)
    assert math.isclose(sum(recorded_dts), 0.3, rel_tol=1e-9)


class _FakeTTY(io.StringIO):
    """io.StringIO reports isatty() == False; ProgressBar's in-place '\\r'
    redraws are gated on isatty(), so exercising that path needs a fake
    terminal that overrides it."""

    def isatty(self) -> bool:
        return True


def test_solve_transient_verbose_does_not_print_a_table_per_timestep(capsys):
    # The whole point: verbose=True on a multi-step transient must not
    # print a per-timestep Newton iteration table (that's the scrolling
    # behavior being replaced) — at most once, for establishing the t=0
    # equilibrium when initial=None.
    network, shaft = _build_dynamic_turboshaft(N0=60000.0)

    history = network.solve_transient(
        duration=0.5, dt=0.05, tol=1e-8, max_iter=400, damping=0.3, verbose=True,
    )

    out = capsys.readouterr().out
    assert out.count("Newton-Raphson solve") == 1  # only the t=0 equilibrium solve
    assert len(history.times) == 11  # 10 steps + t=0, confirms the loop actually ran


def test_solve_transient_verbose_with_explicit_initial_prints_no_iteration_table(capsys):
    # With initial= given, solve_transient() never calls Network.solve()
    # outside the (always-quiet) per-timestep loop, so there should be zero
    # "Newton-Raphson solve" lines even though several steps ran.
    network, shaft = _build_dynamic_turboshaft(N0=60000.0)
    initial = _off_equilibrium_initial(network, shaft, N0_target=55000.0)

    network.solve_transient(
        duration=0.2, dt=0.05, initial=initial, tol=1e-8, max_iter=400, damping=0.3,
        verbose=True,
    )

    out = capsys.readouterr().out
    assert "Newton-Raphson solve" not in out


def test_solve_transient_verbose_prints_done_summary_on_completion(capsys):
    network, shaft = _build_dynamic_turboshaft(N0=60000.0)

    history = network.solve_transient(
        duration=0.3, dt=0.1, tol=1e-8, max_iter=400, damping=0.3, verbose=True,
    )

    out = capsys.readouterr().out
    assert "Done:" in out
    assert "3 steps" in out
    assert math.isclose(history.times[-1], 0.3, rel_tol=1e-9)


def test_solve_transient_quiet_by_default_prints_nothing(capsys):
    network, shaft = _build_dynamic_turboshaft(N0=60000.0)

    network.solve_transient(duration=0.2, dt=0.05, tol=1e-8, max_iter=400, damping=0.3)

    assert capsys.readouterr().out == ""


def test_solve_transient_verbose_interactive_redraws_progress_bar_in_place(monkeypatch):
    # On a real terminal, verbose=True should redraw one fixed line via '\r'
    # for every accepted step (never a newline until the run finishes), and
    # the final line should be colored green. initial= is given explicitly
    # so the only verbose output is the transient bar's own — no Newton bar
    # from establishing a t=0 equilibrium mixed in.
    network, shaft = _build_dynamic_turboshaft(N0=60000.0)
    initial = _off_equilibrium_initial(network, shaft, N0_target=55000.0)
    fake_stdout = _FakeTTY()
    monkeypatch.setattr("sys.stdout", fake_stdout)

    network.solve_transient(
        duration=0.3, dt=0.1, initial=initial, tol=1e-8, max_iter=400, damping=0.3,
        verbose=True,
    )

    out = fake_stdout.getvalue()
    # One redraw per accepted step (3) plus the final colored finish — all
    # '\r'-prefixed, never a '\n' until the very end.
    assert out.count("\r") == 4  # 3 in-loop renders + 1 finish
    assert out.count("\n") == 1  # only the final newline
    assert "\033[32m" in out  # green on completion
    assert "Done:" in out
    assert "step 3" in out
