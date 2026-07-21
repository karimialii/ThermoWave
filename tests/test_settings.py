import math

from thermowave.core.settings import settings


def test_default_units_are_si():
    assert settings.pressure_unit == "Pa"
    assert settings.temperature_unit == "K"


def test_pressure_to_si_default_pa_is_passthrough():
    assert settings.pressure_to_si(101325.0) == 101325.0


def test_pressure_to_si_converts_bar():
    settings.pressure_unit = "bar"
    assert math.isclose(settings.pressure_to_si(1.013), 101300.0, rel_tol=1e-9)


def test_pressure_to_si_converts_atm():
    settings.pressure_unit = "atm"
    assert math.isclose(settings.pressure_to_si(1.0), 101325.0, rel_tol=1e-9)


def test_temperature_to_si_default_kelvin_is_passthrough():
    assert settings.temperature_to_si(300.0) == 300.0


def test_temperature_to_si_converts_celsius():
    settings.temperature_unit = "C"
    assert math.isclose(settings.temperature_to_si(15.0), 288.15, rel_tol=1e-9)
