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

supports_two_phase()/supports_entropy() are typed as TypeIs (PEP 742, via
typing_extensions since the package's own floor is Python 3.10) so mypy can
actually narrow `fluid: BaseFluid` to the matching Protocol below -- but
ONLY inside an `if`/`assert` expression that calls one of them directly, in
the same function scope; a bare statement (the require_*() call below on its
own line) gives mypy no narrowing at all, since it doesn't look inside a
-> None helper to see it would have raised. Every call site in this
codebase that calls require_two_phase()/require_entropy() as a bare
statement therefore follows it with a redundant `assert supports_*(fluid)`
-- unreachable unless require_*() already raised, so stripping `assert`
under `python -O` changes nothing observable (the real, unconditional
enforcement stays in require_*()'s own raise); the assert exists purely to
give mypy the narrowing signal `require_*()` itself can't provide.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

from typing_extensions import TypeIs

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


class TwoPhaseFluid(Protocol):
    """Structural type matching _TWO_PHASE_METHODS, PLUS BaseFluid's own 5
    abstract methods (every real two-phase fluid implements BaseFluid fully
    -- redeclared here, not inherited, since Protocol structural matching
    doesn't need or want a real base-class relationship) so a function
    typed against this Protocol can call density_ph/temperature_ph/etc.
    too, not just the two-phase-specific methods -- see this module's own
    docstring for why this is a Protocol (checked structurally) rather than
    a real base class every two-phase-capable fluid inherits from."""

    name: str

    def density_ph(self, P: float, h: float) -> float: ...
    def temperature_ph(self, P: float, h: float) -> float: ...
    def enthalpy_pt(self, P: float, T: float) -> float: ...
    def cp(self, P: float, T: float) -> float: ...
    def cv(self, P: float, T: float) -> float: ...

    def saturation_temperature(self, P: float) -> float: ...
    def saturated_liquid_enthalpy(self, P: float) -> float: ...
    def saturated_vapor_enthalpy(self, P: float) -> float: ...
    def enthalpy_pq(self, P: float, x: float) -> float: ...
    def quality_ph(self, P: float, h: float) -> float: ...


class EntropyFluid(Protocol):
    """Structural type matching _ENTROPY_METHODS."""

    name: str

    def entropy_ph(self, P: float, h: float) -> float: ...
    def enthalpy_ps(self, P: float, s: float) -> float: ...


def supports_two_phase(fluid: "BaseFluid") -> TypeIs[TwoPhaseFluid]:
    return all(hasattr(fluid, m) for m in _TWO_PHASE_METHODS)


def supports_entropy(fluid: "BaseFluid") -> TypeIs[EntropyFluid]:
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
