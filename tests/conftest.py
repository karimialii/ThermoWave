import pytest

from thermowave.core.settings import settings


@pytest.fixture(autouse=True)
def reset_settings():
    """Ensure the module-level settings singleton doesn't leak between tests."""
    original_pressure_unit = settings.pressure_unit
    original_temperature_unit = settings.temperature_unit
    yield
    settings.pressure_unit = original_pressure_unit
    settings.temperature_unit = original_temperature_unit
