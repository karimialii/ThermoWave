"""Duck-typed capability checks for two-phase / entropy-aware fluids.

Saturation and quality methods live only on CoolPropFluid (the only fluid
model here that can physically represent phase change) -- they are
deliberately NOT on the BaseFluid abstract interface, since
IdealGasFluid/CanteraFluid have no saturation dome to compute them from.
Entropy methods (entropy_ph/enthalpy_ps) are different: CoolPropFluid has
them, and so does every ConstantCpFluid subclass (IdealGasFluid,
IdealGasMixtureFluid) via a closed-form ideal-gas relation (see
ConstantCpFluid.entropy_ph's docstring) -- only CanteraFluid still lacks
them. Components that need either capability check for it structurally
(hasattr on the specific method names below) rather than with
isinstance(fluid, CoolPropFluid): checking for the actual methods a
component calls, rather than a specific class, means any fluid exposing
that contract qualifies -- including a future REFPROP-backed fluid, or any
other model that computes saturation/entropy properties some other way.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from thermowave.fluids.base_fluid import BaseFluid

# Methods a fluid must expose to back an evaporator/condenser/drum.
_TWO_PHASE_METHODS = (
    "saturation_temperature",
    "saturated_liquid_enthalpy",
    "saturated_vapor_enthalpy",
    "enthalpy_pq",
    "quality_ph",
)

# Methods a fluid must expose for an isentropic-path component (Pump, SteamTurbine).
_ENTROPY_METHODS = ("entropy_ph", "enthalpy_ps")


def supports_two_phase(fluid: "BaseFluid") -> bool:
    return all(hasattr(fluid, m) for m in _TWO_PHASE_METHODS)


def supports_entropy(fluid: "BaseFluid") -> bool:
    return all(hasattr(fluid, m) for m in _ENTROPY_METHODS)


def require_two_phase(fluid: "BaseFluid", component_name: str) -> None:
    if not supports_two_phase(fluid):
        raise ValueError(
            f"{component_name!r} needs a two-phase-capable fluid exposing "
            f"{_TWO_PHASE_METHODS} (e.g. CoolPropFluid) -- got "
            f"{type(fluid).__name__} {getattr(fluid, 'name', '?')!r}, which has no "
            f"saturation/quality model. IdealGasFluid/CanteraFluid cannot boil or "
            f"condense; use CoolPropFluid for phase-change components."
        )


def require_entropy(fluid: "BaseFluid", component_name: str) -> None:
    if not supports_entropy(fluid):
        raise ValueError(
            f"{component_name!r} needs a fluid exposing {_ENTROPY_METHODS} for its "
            f"isentropic path (CoolPropFluid, IdealGasFluid, or IdealGasMixtureFluid "
            f"all have it) -- got {type(fluid).__name__} {getattr(fluid, 'name', '?')!r}. "
            f"The Cantera model doesn't provide entropy here."
        )
