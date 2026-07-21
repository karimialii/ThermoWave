import math

import pytest

from thermowave.components.generator import Generator
from thermowave.components.simple_generator import SimpleGenerator
from thermowave.maps.torque_speed_map import TorqueSpeedMap

_MAP_PATH = "tests/fixtures/simple_generator_map.gen"


class _FakeShaftComponent:
    """Stand-in for a Turbine/Compressor: report_metrics() only."""

    name = "shaft1"

    def __init__(self, metrics):
        self._metrics = metrics

    def report_metrics(self, state):
        return self._metrics


# --- TorqueSpeedMap -----------------------------------------------------


def test_torque_speed_map_from_text_parses_pairs_and_skips_comments():
    text = "# header\n40000 20.0\n\n60000 24.0\n"
    tsm = TorqueSpeedMap.from_text(text)
    assert math.isclose(tsm.torque(40000.0), 20.0)
    assert math.isclose(tsm.torque(60000.0), 24.0)


def test_torque_speed_map_interpolates_between_points():
    tsm = TorqueSpeedMap.from_text("40000 20.0\n60000 24.0\n")
    assert math.isclose(tsm.torque(50000.0), 22.0, rel_tol=1e-9)


def test_torque_speed_map_clamps_outside_range():
    tsm = TorqueSpeedMap.from_text("40000 20.0\n60000 24.0\n")
    assert math.isclose(tsm.torque(10000.0), 20.0)
    assert math.isclose(tsm.torque(100000.0), 24.0)


def test_torque_speed_map_mid_speed():
    tsm = TorqueSpeedMap.from_text("40000 20.0\n60000 24.0\n")
    assert math.isclose(tsm.mid_speed(), 50000.0)


def test_torque_speed_map_from_file_matches_fixture():
    tsm = TorqueSpeedMap.from_file(_MAP_PATH)
    assert math.isclose(tsm.torque(50000.0), 22.0, rel_tol=1e-9)


# --- SimpleGenerator ------------------------------------------------------


def test_simple_generator_ports_and_category():
    shaft = _FakeShaftComponent({"power [W]": 1000.0})
    gen = SimpleGenerator(name="g1", component=shaft, efficiency=0.9)
    assert gen.ports() == {}
    assert gen.report_category() == "generator"
    assert gen.residuals(state=None) == []


def test_simple_generator_scales_shaft_power_by_efficiency():
    shaft = _FakeShaftComponent({"power [W]": 1000.0, "N [rev/min]": 50000.0})
    gen = SimpleGenerator(name="g1", component=shaft, efficiency=0.9)
    metrics = gen.report_metrics(state=None)
    assert math.isclose(metrics["power [W]"], 900.0)
    assert math.isclose(metrics["eta [-]"], 0.9)
    assert math.isclose(metrics["N [rev/min]"], 50000.0)


def test_simple_generator_raises_if_shaft_reports_no_power():
    shaft = _FakeShaftComponent({"eta_s [-]": 0.8})
    gen = SimpleGenerator(name="g1", component=shaft, efficiency=0.9)
    with pytest.raises(ValueError, match="power \\[W\\]"):
        gen.report_metrics(state=None)


# --- Generator (map-based) -------------------------------------------------


def test_generator_ports_and_category():
    shaft = _FakeShaftComponent({"N [rev/min]": 50000.0})
    gen = Generator(name="g1", component=shaft, map_path=_MAP_PATH)
    assert gen.ports() == {}
    assert gen.report_category() == "generator"
    assert gen.residuals(state=None) == []


def test_generator_power_matches_torque_times_omega_hand_calc():
    shaft = _FakeShaftComponent({"N [rev/min]": 50000.0})
    gen = Generator(name="g1", component=shaft, map_path=_MAP_PATH)
    torque = 22.0  # from the fixture map at 50000 rev/min
    omega = 50000.0 * 2.0 * math.pi / 60.0
    expected_power = torque * omega
    metrics = gen.report_metrics(state=None)
    assert math.isclose(metrics["power [W]"], expected_power, rel_tol=1e-9)
    assert math.isclose(metrics["N [rev/min]"], 50000.0)
    assert math.isclose(metrics["eta [-]"], 1.0)


def test_generator_applies_efficiency_on_top_of_map():
    shaft = _FakeShaftComponent({"N [rev/min]": 50000.0})
    gen = Generator(name="g1", component=shaft, map_path=_MAP_PATH, efficiency=0.95)
    metrics_full = Generator(name="g0", component=shaft, map_path=_MAP_PATH).report_metrics(None)
    metrics = gen.report_metrics(state=None)
    assert math.isclose(metrics["power [W]"], metrics_full["power [W]"] * 0.95, rel_tol=1e-9)


def test_generator_raises_if_shaft_reports_no_speed():
    shaft = _FakeShaftComponent({"power [W]": 1000.0})
    gen = Generator(name="g1", component=shaft, map_path=_MAP_PATH)
    with pytest.raises(ValueError, match="N \\[rev/min\\]"):
        gen.report_metrics(state=None)
