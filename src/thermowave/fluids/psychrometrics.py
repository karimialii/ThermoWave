"""Duck-typed capability check for humid-air-aware fluids.

Relative-humidity/wet-bulb/dew-point methods live only on HumidAirFluid (the
only fluid model here with a humidity concept at all) -- they are
deliberately NOT on the BaseFluid abstract interface, since no other fluid
model here (IdealGasFluid, CoolPropFluid, CanteraFluid, ...) has a humidity
ratio to compute them from. Components that need this capability check for
it structurally (hasattr on the specific method names below) rather than
with isinstance(fluid, HumidAirFluid): checking for the actual methods a
component calls, rather than a specific class, means any fluid exposing
that contract qualifies -- same rationale as thermowave.fluids.two_phase.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from thermowave.fluids.base_fluid import BaseFluid

# Methods a fluid must expose to back a psychrometric/HVAC component.
_HUMID_AIR_METHODS = (
    "relative_humidity_pt",
    "wet_bulb_pt",
    "dew_point_pt",
)


def supports_humid_air(fluid: "BaseFluid") -> bool:
    return all(hasattr(fluid, m) for m in _HUMID_AIR_METHODS)


def require_humid_air(fluid: "BaseFluid", component_name: str) -> None:
    if not supports_humid_air(fluid):
        raise ValueError(
            f"{component_name!r} needs a humid-air-capable fluid exposing "
            f"{_HUMID_AIR_METHODS} (e.g. HumidAirFluid) -- got "
            f"{type(fluid).__name__} {getattr(fluid, 'name', '?')!r}, which has no "
            f"psychrometric model."
        )
