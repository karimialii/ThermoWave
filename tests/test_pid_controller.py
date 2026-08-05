import math

import pytest

from thermowave.components.compressor import Compressor
from thermowave.components.pid_controller import PIDController
from thermowave.components.sensor import Sensor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)
GAMMA = 1005.0 / (1005.0 - 287.05)


class _FakeState:
    def __init__(self, fluid, node_values, params=None, node_N=None):
        self.fluid = fluid
        self._node_values = node_values
        self._params = params or {}
        self._node_N = node_N or {}

    def node(self, name):
        return self._node_values[name]

    def fluid_at(self, name):
        return self.fluid

    def mdot(self, name):
        raise KeyError(name)

    def param(self, name):
        return self._params[name]

    def N(self, name):
        return self._node_N[name]


def _make_pid(**overrides):
    comp = Compressor(name="c1", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None)
    sensor = Sensor(name="s1")
    kwargs = dict(
        name="pid1", sensor=sensor, quantity="T [K]", component=comp, free_param="shaft",
        setpoint=420.0, Kp=60.0, Ki=50.0, Kd=0.0, output0=60000.0,
    )
    kwargs.update(overrides)
    return PIDController(**kwargs), comp, sensor


def _temp_state(fluid, node, P, T, params=None):
    h = fluid.enthalpy_pt(P, T)
    return _FakeState(fluid=fluid, node_values={node: (P, h)}, params=params)


def test_pid_raises_if_free_param_not_declared_free():
    comp = Compressor(name="c1", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=50000.0)
    sensor = Sensor(name="s1")
    with pytest.raises(ValueError, match="doesn't currently declare"):
        PIDController(
            name="pid1", sensor=sensor, quantity="T [K]", component=comp, free_param="shaft",
            setpoint=420.0, Kp=1.0,
        )


def test_pid_ports_is_empty():
    pid, _, _ = _make_pid()
    assert pid.ports() == {}


def test_pid_residual_pins_actual_free_param_to_current_output():
    pid, comp, _ = _make_pid()
    state = _FakeState(fluid=AIR, node_values={}, node_N={"c1.shaft": 55000.0})
    pid.output = 55000.0
    assert math.isclose(pid.residuals(state)[0], 0.0, abs_tol=1e-9)

    pid.output = 55500.0
    assert math.isclose(pid.residuals(state)[0], -500.0, abs_tol=1e-9)


def test_pid_step_proportional_only_moves_output_toward_setpoint():
    pid, comp, sensor = _make_pid(Ki=0.0, Kd=0.0)
    state = _temp_state(AIR, "s1.tap", P=300000.0, T=445.0)  # error = 420 - 445 = -25
    new_output = pid.step(state, dt=0.05)
    # output = bias + Kp*error = 60000 + 60*(-25) = 58500
    assert math.isclose(new_output, 58500.0, rel_tol=1e-9)
    assert pid.output == new_output


def test_pid_step_integral_accumulates_over_successive_calls():
    pid, comp, sensor = _make_pid(Kp=0.0, Ki=100.0, Kd=0.0)
    state = _temp_state(AIR, "s1.tap", P=300000.0, T=430.0)  # error = -10, constant
    pid.step(state, dt=0.1)   # integral = -1.0  -> output = 60000 + 100*-1.0 = 59900
    pid.step(state, dt=0.1)   # integral = -2.0  -> output = 60000 + 100*-2.0 = 59800
    assert math.isclose(pid.output, 59800.0, rel_tol=1e-9)


def test_pid_step_derivative_uses_error_change_over_dt():
    pid, comp, sensor = _make_pid(Kp=0.0, Ki=0.0, Kd=2.0)
    state1 = _temp_state(AIR, "s1.tap", P=300000.0, T=420.0)  # error = 0
    state2 = _temp_state(AIR, "s1.tap", P=300000.0, T=440.0)  # error = -20
    pid.step(state1, dt=0.1)
    pid.step(state2, dt=0.1)
    # derivative = (-20 - 0) / 0.1 = -200 -> output = 60000 + 2*(-200) = 59600
    assert math.isclose(pid.output, 59600.0, rel_tol=1e-9)


def test_pid_step_derivative_ignores_a_setpoint_step_with_constant_measurement():
    """Derivative-on-measurement (not derivative-on-error): a step change in
    self.setpoint with the *measurement* unchanged between calls must not
    inject a derivative kick into output, since the plant hasn't actually
    moved -- only what we're asking of it has.
    """
    pid, comp, sensor = _make_pid(Kp=0.0, Ki=0.0, Kd=5.0)
    state = _temp_state(AIR, "s1.tap", P=300000.0, T=420.0)  # measured = 420, constant
    pid.step(state, dt=0.1)
    pid.setpoint = 400.0  # a step change in the command, not the plant
    new_output = pid.step(state, dt=0.1)
    assert math.isclose(new_output, pid.bias, abs_tol=1e-9)


def test_pid_step_clamps_output_and_stops_integral_windup_when_saturated():
    pid, comp, sensor = _make_pid(Kp=0.0, Ki=1000.0, Kd=0.0, output_min=59000.0, output_max=61000.0)
    state = _temp_state(AIR, "s1.tap", P=300000.0, T=430.0)  # error = -10
    pid.step(state, dt=1.0)  # unclamped would be 60000 - 10000 = 50000, clamps to 59000
    assert pid.output == 59000.0
    # Integral shouldn't have accumulated (discarded while saturated) -> a
    # second identical step gives the same clamped output, not a further drop.
    pid.step(state, dt=1.0)
    assert pid.output == 59000.0


def test_pid_first_step_has_no_derivative_term():
    pid, comp, sensor = _make_pid(Kp=0.0, Ki=0.0, Kd=5.0)
    state = _temp_state(AIR, "s1.tap", P=300000.0, T=440.0)  # error = -20
    pid.step(state, dt=0.1)
    assert math.isclose(pid.output, pid.bias, abs_tol=1e-9)


def test_pid_report_metrics_reflects_setpoint_measured_error_and_output():
    pid, comp, sensor = _make_pid(Kp=0.0, Ki=0.0, Kd=0.0)
    state = _temp_state(AIR, "s1.tap", P=300000.0, T=430.0)
    metrics = pid.report_metrics(state)
    assert math.isclose(metrics["target [-]"], 420.0)
    assert math.isclose(metrics["measured [-]"], 430.0)
    assert math.isclose(metrics["error [-]"], -10.0)
    assert math.isclose(metrics["output [-]"], pid.output)


def test_pid_controller_drives_compressor_outlet_temperature_over_transient():
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    comp = Compressor(name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None)
    sensor = Sensor(name="s1")
    snk = Sink(name="snk")

    target_T = 420.0
    pid = PIDController(
        name="pid", sensor=sensor, quantity="T [K]", component=comp, free_param="shaft",
        setpoint=target_T, Kp=60.0, Ki=50.0, Kd=0.0,
        output0=60000.0, output_min=10000.0, output_max=100000.0,
    )

    network = Network(fluid=AIR)
    for component in (src, comp, sensor, pid, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", sensor, "tap")
    network.connect(comp, "out", snk, "in")

    history = network.solve_transient(
        duration=8.0, dt=0.1, tol=1e-6, max_iter=300, damping=0.5,
    )

    def _T(step):
        return AIR.temperature_ph(step.node_P["comp.out"], step.node_h["comp.out"])

    T_start = _T(history.steps[0])
    T_end = _T(history.steps[-1])
    assert T_start > target_T  # starts hot (N0 too high for the target)
    assert abs(T_end - target_T) < abs(T_start - target_T)  # error shrank
    assert abs(T_end - target_T) < 2.0  # nearly settled after 8 s


def test_pid_feedforward_replaces_the_constant_bias_each_step():
    # The bias a PID trims around can move with the plant instead of being
    # pinned to whatever it was at t=0 -- e.g. fuel biased from commanded
    # load, so the integral term never has to hunt across the whole
    # actuator range to find the operating point at all.
    moving = {"value": 60000.0}
    pid, comp, sensor = _make_pid(
        Kp=0.0, Ki=0.0, Kd=0.0, feedforward=lambda state: moving["value"]
    )
    state = _temp_state(AIR, "s1.tap", P=300000.0, T=420.0)  # error = 0

    pid.step(state, dt=1.0)
    assert math.isclose(pid.output, 60000.0, rel_tol=1e-9)

    moving["value"] = 55000.0
    pid.step(state, dt=1.0)
    assert math.isclose(pid.output, 55000.0, rel_tol=1e-9)


def test_pid_feedforward_is_added_to_the_pid_terms_not_instead_of_them():
    pid, comp, sensor = _make_pid(Kp=2.0, Ki=0.0, Kd=0.0, feedforward=lambda state: 50000.0)
    state = _temp_state(AIR, "s1.tap", P=300000.0, T=430.0)  # error = -10

    pid.step(state, dt=1.0)

    assert math.isclose(pid.output, 50000.0 + 2.0 * -10.0, rel_tol=1e-9)


def test_pid_bias_is_public_so_a_schedule_can_drive_it():
    pid, comp, sensor = _make_pid(Kp=0.0, Ki=0.0, Kd=0.0, output0=60000.0)
    state = _temp_state(AIR, "s1.tap", P=300000.0, T=420.0)  # error = 0

    pid.bias = 58000.0
    pid.step(state, dt=1.0)

    assert math.isclose(pid.output, 58000.0, rel_tol=1e-9)


def test_pid_reset_clears_integral_windup_and_previous_error():
    # A PIDController carries this state on the instance, so reusing one
    # across two solve_transient() runs would otherwise start the second
    # with the first's wound-up integral.
    pid, comp, sensor = _make_pid(Kp=0.0, Ki=1000.0, Kd=0.0, output0=60000.0)
    state = _temp_state(AIR, "s1.tap", P=300000.0, T=430.0)  # error = -10
    pid.step(state, dt=1.0)
    assert not math.isclose(pid.output, 60000.0, rel_tol=1e-9)  # integral moved it

    pid.reset()

    assert pid.output == 60000.0
    assert pid.bias == 60000.0
    assert pid._integral == 0.0
    assert pid._prev_measured is None
    # And a fresh step behaves like the very first one again.
    pid.step(state, dt=1.0)
    assert math.isclose(pid.output, 60000.0 + 1000.0 * -10.0, rel_tol=1e-9)


def test_pid_reset_can_reseed_the_operating_point():
    pid, comp, sensor = _make_pid(Kp=0.0, Ki=0.0, Kd=0.0, output0=60000.0)

    pid.reset(output0=52000.0)

    assert pid.output == 52000.0
    assert pid.bias == 52000.0
