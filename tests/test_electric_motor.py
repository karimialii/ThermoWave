import math

import pytest

from thermowave.components.electric_motor import ElectricMotor


class _FakeMechComponent:
    """Stand-in for a Compressor/pump: report_metrics() only."""

    name = "comp1"

    def __init__(self, metrics):
        self._metrics = metrics

    def report_metrics(self, state):
        return self._metrics


def test_electric_motor_rejects_efficiency_out_of_range():
    comp = _FakeMechComponent({"power [W]": 1000.0})
    with pytest.raises(ValueError, match="efficiency"):
        ElectricMotor(name="m1", component=comp, efficiency=0.0)
    with pytest.raises(ValueError, match="efficiency"):
        ElectricMotor(name="m1", component=comp, efficiency=1.5)


def test_electric_motor_ports_and_category():
    comp = _FakeMechComponent({"power [W]": 1000.0})
    motor = ElectricMotor(name="m1", component=comp, efficiency=0.9)
    assert motor.ports() == {}
    assert motor.report_category() == "motor"
    assert motor.residuals(state=None) == []


def test_electric_motor_divides_required_shaft_power_by_efficiency():
    comp = _FakeMechComponent({"power [W]": 900.0, "N [rev/min]": 40000.0})
    motor = ElectricMotor(name="m1", component=comp, efficiency=0.9)
    metrics = motor.report_metrics(state=None)
    assert math.isclose(metrics["power [W]"], 1000.0)
    assert math.isclose(metrics["eta [-]"], 0.9)
    assert math.isclose(metrics["N [rev/min]"], 40000.0)


def test_electric_motor_omits_speed_when_mechanical_component_has_none():
    comp = _FakeMechComponent({"power [W]": 900.0})
    motor = ElectricMotor(name="m1", component=comp, efficiency=0.9)
    metrics = motor.report_metrics(state=None)
    assert "N [rev/min]" not in metrics


def test_electric_motor_raises_if_mechanical_component_reports_no_power():
    comp = _FakeMechComponent({"eta_s [-]": 0.8})
    motor = ElectricMotor(name="m1", component=comp, efficiency=0.9)
    with pytest.raises(ValueError, match="power \\[W\\]"):
        motor.report_metrics(state=None)


def test_electric_motor_formula_is_the_algebraic_inverse_of_simple_generator():
    # SimpleGenerator: power_elec = shaft_power * eta.
    # ElectricMotor:    power_elec = shaft_power_required / eta.
    # Feeding one's output back in as the other's "required power" input at
    # the same eta should recover the original value exactly — confirms the
    # two formulas are genuine inverses of each other, not just similarly
    # shaped.
    from thermowave.components.simple_generator import SimpleGenerator

    shaft_power = 1000.0
    eta = 0.9
    comp = _FakeMechComponent({"power [W]": shaft_power})
    gen = SimpleGenerator(name="g1", component=comp, efficiency=eta)
    elec_out = gen.report_metrics(state=None)["power [W]"]

    driven = _FakeMechComponent({"power [W]": elec_out})
    motor = ElectricMotor(name="m1", component=driven, efficiency=eta)
    round_trip = motor.report_metrics(state=None)["power [W]"]
    assert math.isclose(round_trip, shaft_power, rel_tol=1e-9)
