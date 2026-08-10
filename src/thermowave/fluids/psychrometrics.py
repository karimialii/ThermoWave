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

supports_humid_air() is TypeIs-typed for the same reason and with the same
call-site caveat as thermowave.fluids.two_phase.supports_two_phase() -- see
that module's docstring; a bare require_humid_air() statement gives mypy no
narrowing, so callers follow it with a redundant `assert supports_humid_air(fluid)`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

from typing_extensions import TypeIs

if TYPE_CHECKING:
    from thermowave.fluids.base_fluid import BaseFluid

# Methods a fluid must expose to back a psychrometric/HVAC component.
_HUMID_AIR_METHODS = (
    "relative_humidity_pt",
    "wet_bulb_pt",
    "dew_point_pt",
)


class HumidAirCapableFluid(Protocol):
    """Structural type matching _HUMID_AIR_METHODS."""

    name: str
    W: float

    def relative_humidity_pt(self, P: float, T: float) -> float: ...
    def wet_bulb_pt(self, P: float, T: float) -> float: ...
    def dew_point_pt(self, P: float, T: float) -> float: ...


def supports_humid_air(fluid: "BaseFluid") -> TypeIs[HumidAirCapableFluid]:
    return all(hasattr(fluid, m) for m in _HUMID_AIR_METHODS)


def require_humid_air(fluid: "BaseFluid", component_name: str) -> None:
    if not supports_humid_air(fluid):
        raise ValueError(
            f"{component_name!r} needs a humid-air-capable fluid exposing "
            f"{_HUMID_AIR_METHODS} (e.g. HumidAirFluid) -- got "
            f"{type(fluid).__name__} {getattr(fluid, 'name', '?')!r}, which has no "
            f"psychrometric model."
        )
