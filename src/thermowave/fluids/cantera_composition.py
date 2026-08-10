"""Duck-typed capability check for Cantera-flavored (composition-blendable)
fluids.

mass_fractions()/mechanism live only on CanteraFluid and
_CanteraCompositionFluid (Combustor's own equilibrium-product fluid) --
deliberately NOT on the BaseFluid abstract interface, since no other fluid
model here has a per-species composition to expose. Junction.merge_fluids()
and Combustor._equilibrium_mixture() both need to know whether an inlet
fluid can be mass-weighted/mixed into a Cantera Solution before doing so;
this is the same structural check both already made independently via
inline hasattr() calls, collected here into one shared module -- same
rationale, and same TypeIs call-site caveat, as
thermowave.fluids.two_phase/thermowave.fluids.psychrometrics: TypeIs only
narrows inside the `if`/`assert` expression that calls it directly, in the
same function scope.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

from typing_extensions import TypeIs

if TYPE_CHECKING:
    from thermowave.fluids.base_fluid import BaseFluid

# Attributes/methods a fluid must expose to be blended/passed as a Cantera
# composition source. "mechanism" is a plain attribute, not a method --
# hasattr() covers both the same way.
_CANTERA_COMPOSITION_ATTRS = ("mass_fractions", "mechanism")


class CanteraCompositionFluid(Protocol):
    """Structural type matching _CANTERA_COMPOSITION_ATTRS."""

    name: str
    mechanism: str

    def mass_fractions(self) -> dict[str, float]: ...


def supports_cantera_composition(fluid: "BaseFluid") -> TypeIs[CanteraCompositionFluid]:
    return all(hasattr(fluid, a) for a in _CANTERA_COMPOSITION_ATTRS)


def require_cantera_composition(fluid: "BaseFluid", component_name: str) -> None:
    if not supports_cantera_composition(fluid):
        raise ValueError(
            f"{component_name!r} needs a Cantera-composition-capable fluid exposing "
            f"{_CANTERA_COMPOSITION_ATTRS} (e.g. CanteraFluid) -- got "
            f"{type(fluid).__name__} {getattr(fluid, 'name', '?')!r}."
        )
