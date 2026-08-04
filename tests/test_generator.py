import math

from thermowave.components.generator import Generator
from thermowave.components.simple_generator import SimpleGenerator
from thermowave.maps.torque_speed_map import TorqueSpeedMap

_MAP_PATH = "tests/fixtures/simple_generator_map.gen"


class _FakeState:
    """Stand-in for core.network.NetworkState -- power/N read via ports."""

    def __init__(self, node_N=None, node_signal=None):
        self._node_N = node_N or {}
        self._node_signal = node_signal or {}

    def N(self, name):
        return self._node_N[name]

    def signal(self, name):
        return self._node_signal[name]


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
    gen = SimpleGenerator(name="g1", efficiency=0.9)
    assert gen.ports() == {}
    assert gen.report_category() == "generator"
    assert gen.residuals(state=None) == []


def test_simple_generator_scales_shaft_power_by_efficiency():
    gen = SimpleGenerator(name="g1", efficiency=0.9)
    state = _FakeState(
        node_N={gen._shaft_port: 50000.0}, node_signal={gen._power_port: 1000.0}
    )
    metrics = gen.report_metrics(state)
    assert math.isclose(metrics["power [W]"], 900.0)
    assert math.isclose(metrics["eta [-]"], 0.9)
    assert math.isclose(metrics["N [rev/min]"], 50000.0)


def test_simple_generator_omits_speed_when_shaft_not_connected():
    gen = SimpleGenerator(name="g1", efficiency=0.9)
    state = _FakeState(node_signal={gen._power_port: 1000.0})
    metrics = gen.report_metrics(state)
    assert "N [rev/min]" not in metrics


# --- Generator (map-based) -------------------------------------------------


def test_generator_ports_and_category():
    gen = Generator(name="g1", map_path=_MAP_PATH)
    assert gen.ports() == {}
    assert gen.report_category() == "generator"
    assert gen.residuals(state=None) == []


def test_generator_power_matches_torque_times_omega_hand_calc():
    gen = Generator(name="g1", map_path=_MAP_PATH)
    state = _FakeState(node_N={gen._shaft_port: 50000.0})
    torque = 22.0  # from the fixture map at 50000 rev/min
    omega = 50000.0 * 2.0 * math.pi / 60.0
    expected_power = torque * omega
    metrics = gen.report_metrics(state)
    assert math.isclose(metrics["power [W]"], expected_power, rel_tol=1e-9)
    assert math.isclose(metrics["N [rev/min]"], 50000.0)
    assert math.isclose(metrics["eta [-]"], 1.0)


def test_generator_applies_efficiency_on_top_of_map():
    gen = Generator(name="g1", map_path=_MAP_PATH, efficiency=0.95)
    gen_full = Generator(name="g0", map_path=_MAP_PATH)
    state = _FakeState(node_N={gen._shaft_port: 50000.0, gen_full._shaft_port: 50000.0})
    metrics_full = gen_full.report_metrics(state)
    metrics = gen.report_metrics(state)
    assert math.isclose(metrics["power [W]"], metrics_full["power [W]"] * 0.95, rel_tol=1e-9)


def test_generator_provided_signal_matches_report_metrics_power():
    gen = Generator(name="g1", map_path=_MAP_PATH)
    state = _FakeState(node_N={gen._shaft_port: 50000.0})
    assert math.isclose(
        gen.provided_signal_values(state)[gen._power_port],
        gen.report_metrics(state)["power [W]"],
    )
