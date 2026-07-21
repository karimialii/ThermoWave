import pytest

from thermowave.components.heat_transfer import ThermalMass


class _FakeState:
    """Minimal stand-in for core.network.NetworkState."""

    def __init__(self, params):
        self._params = params

    def param(self, name: str) -> float:
        return self._params[name]


class _FakePath:
    """Stand-in heat-path component: Q(state) returns a fixed value."""

    def __init__(self, name, Q_value):
        self.name = name
        self._Q_value = Q_value

    def Q(self, state):
        return self._Q_value


def test_thermal_mass_rejects_non_positive_capacitance():
    with pytest.raises(ValueError, match="thermal_capacitance must be > 0"):
        ThermalMass(name="tm1", thermal_capacitance=0.0, T0=300.0)


def test_thermal_mass_ports_and_residuals_are_empty():
    tm = ThermalMass(name="tm1", thermal_capacitance=500.0, T0=300.0)
    assert tm.ports() == {}
    assert tm.residuals(_FakeState({})) == []


def test_thermal_mass_differential_parameters_seeded_from_t0():
    tm = ThermalMass(name="tm1", thermal_capacitance=500.0, T0=310.0)
    assert tm.differential_parameters() == {"T": 310.0}


def test_thermal_mass_state_derivative_zero_with_no_heat_sources():
    tm = ThermalMass(name="tm1", thermal_capacitance=500.0, T0=300.0)
    state = _FakeState({"tm1.T": 300.0})
    assert tm.state_derivative(state) == {"T": 0.0}


def test_thermal_mass_state_derivative_single_source():
    tm = ThermalMass(name="tm1", thermal_capacitance=500.0, T0=300.0)
    path = _FakePath("p1", Q_value=1000.0)
    tm.heat_sources = [(path, 1.0)]
    state = _FakeState({"tm1.T": 300.0})
    derivative = tm.state_derivative(state)
    assert derivative == {"T": 1000.0 / 500.0}


def test_thermal_mass_state_derivative_sums_signed_sources():
    tm = ThermalMass(name="tm1", thermal_capacitance=200.0, T0=300.0)
    gain = _FakePath("gain", Q_value=800.0)
    loss = _FakePath("loss", Q_value=300.0)
    tm.heat_sources = [(gain, 1.0), (loss, -1.0)]
    state = _FakeState({"tm1.T": 300.0})
    derivative = tm.state_derivative(state)
    assert derivative == {"T": (800.0 - 300.0) / 200.0}


def test_thermal_mass_report_metrics_reflects_current_t_and_net_q():
    tm = ThermalMass(name="tm1", thermal_capacitance=500.0, T0=300.0)
    path = _FakePath("p1", Q_value=1000.0)
    tm.heat_sources = [(path, 1.0)]
    state = _FakeState({"tm1.T": 345.0})
    metrics = tm.report_metrics(state)
    assert metrics["T [K]"] == 345.0
    assert metrics["Q_net [W]"] == 1000.0


def test_thermal_mass_report_category_is_thermal_mass():
    tm = ThermalMass(name="tm1", thermal_capacitance=500.0, T0=300.0)
    assert tm.report_category() == "thermal_mass"
