import math

import pytest

from thermowave.components.electric_motor import ElectricMotor
from thermowave.components.simple_generator import SimpleGenerator


class _FakeState:
    def __init__(self, node_N=None, node_signal=None):
        self._node_N = node_N or {}
        self._node_signal = node_signal or {}

    def N(self, name):
        return self._node_N[name]

    def signal(self, name):
        return self._node_signal[name]


def test_electric_motor_rejects_efficiency_out_of_range():
    with pytest.raises(ValueError, match="efficiency"):
        ElectricMotor(name="m1", efficiency=0.0)
    with pytest.raises(ValueError, match="efficiency"):
        ElectricMotor(name="m1", efficiency=1.5)


def test_electric_motor_ports_and_category():
    motor = ElectricMotor(name="m1", efficiency=0.9)
    assert motor.ports() == {}
    assert motor.report_category() == "motor"
    assert motor.residuals(state=None) == []


def test_electric_motor_divides_required_shaft_power_by_efficiency():
    motor = ElectricMotor(name="m1", efficiency=0.9)
    state = _FakeState(
        node_N={motor._shaft_port: 40000.0}, node_signal={motor._power_port: 900.0}
    )
    metrics = motor.report_metrics(state)
    assert math.isclose(metrics["power [W]"], 1000.0)
    assert math.isclose(metrics["eta [-]"], 0.9)
    assert math.isclose(metrics["N [rev/min]"], 40000.0)


def test_electric_motor_omits_speed_when_shaft_not_connected():
    motor = ElectricMotor(name="m1", efficiency=0.9)
    state = _FakeState(node_signal={motor._power_port: 900.0})
    metrics = motor.report_metrics(state)
    assert "N [rev/min]" not in metrics


def test_electric_motor_formula_is_the_algebraic_inverse_of_simple_generator():
    # SimpleGenerator: power_elec = shaft_power * eta.
    # ElectricMotor:    power_elec = shaft_power_required / eta.
    # Feeding one's output back in as the other's "required power" input at
    # the same eta should recover the original value exactly — confirms the
    # two formulas are genuine inverses of each other, not just similarly
    # shaped.
    shaft_power = 1000.0
    eta = 0.9
    gen = SimpleGenerator(name="g1", efficiency=eta)
    gen_state = _FakeState(node_signal={gen._power_port: shaft_power})
    elec_out = gen.report_metrics(gen_state)["power [W]"]

    motor = ElectricMotor(name="m1", efficiency=eta)
    motor_state = _FakeState(node_signal={motor._power_port: elec_out})
    round_trip = motor.report_metrics(motor_state)["power [W]"]
    assert math.isclose(round_trip, shaft_power, rel_tol=1e-9)
