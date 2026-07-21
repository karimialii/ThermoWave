from thermowave.core.constants import PA_PER_BAR, STANDARD_ATMOSPHERE_PA

_PRESSURE_TO_PA = {
    "Pa": 1.0,
    "kPa": 1.0e3,
    "MPa": 1.0e6,
    "bar": PA_PER_BAR,
    "atm": STANDARD_ATMOSPHERE_PA,
}


class Settings:
    """Configurable I/O units for constructing components and displaying results.

    All internal solver/component math is strict SI regardless of these settings;
    conversion happens only at construction time (see Source.__init__).
    """

    def __init__(self) -> None:
        self.pressure_unit = "Pa"
        self.temperature_unit = "K"

    def pressure_to_si(self, value: float) -> float:
        if self.pressure_unit not in _PRESSURE_TO_PA:
            raise ValueError(f"Unsupported pressure_unit: {self.pressure_unit!r}")
        return value * _PRESSURE_TO_PA[self.pressure_unit]

    def temperature_to_si(self, value: float) -> float:
        # Only K/C are supported, unlike pressure's five units above — this
        # is deliberate, not an oversight: every fluid model and component
        # in this codebase works in fixed-point Kelvin only where an offset
        # matters (e.g. equilibrium chemistry, ideal-gas relations), so
        # Fahrenheit's non-trivial affine conversion (and mixed-unit
        # confusion risk) hasn't been worth adding without a concrete use
        # case. Add an "F" branch here (value - 32) * 5/9 + 273.15 if one
        # comes up.
        if self.temperature_unit == "K":
            return value
        if self.temperature_unit == "C":
            return value + 273.15
        raise ValueError(f"Unsupported temperature_unit: {self.temperature_unit!r}")


settings = Settings()
