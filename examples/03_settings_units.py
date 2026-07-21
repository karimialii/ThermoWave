"""settings: module-level singleton controlling non-SI input units.

All internal solver math stays strict SI; `settings` only affects how values
passed into components like Source are converted at construction time.

Run: .venv/bin/python examples/03_settings_units.py
"""

from thermowave.core.settings import settings

print(f"default pressure_unit    : {settings.pressure_unit!r}")
print(f"default temperature_unit : {settings.temperature_unit!r}")
print(f"101325 Pa -> SI          : {settings.pressure_to_si(101325.0):.1f} Pa")
print()

settings.pressure_unit = "bar"
settings.temperature_unit = "C"
print(f"pressure_unit = 'bar', 1.013 bar -> SI : {settings.pressure_to_si(1.013):.1f} Pa")
print(f"temperature_unit = 'C', 26.85 C -> SI  : {settings.temperature_to_si(26.85):.2f} K")

settings.pressure_unit = "atm"
print(f"pressure_unit = 'atm', 1.0 atm -> SI   : {settings.pressure_to_si(1.0):.1f} Pa")

# Reset to SI defaults so other code (or a later example run) isn't affected.
settings.pressure_unit = "Pa"
settings.temperature_unit = "K"
