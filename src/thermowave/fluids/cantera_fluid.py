from __future__ import annotations

from thermowave.core.constants import STANDARD_ATMOSPHERE_PA
from thermowave.fluids.base_fluid import BaseFluid

_REFERENCE_T = 300.0  # K, only used as a starting point before solving HP for T


class _CanteraCompositionFluid(BaseFluid):
    """Lightweight BaseFluid view over one already-known, fixed mass-fraction
    composition, backed by a Cantera Solution the *caller* already loaded
    once (not created here) -- unlike CanteraFluid, this does not parse a
    mechanism file itself, so constructing one of these (once per Newton
    residual evaluation, which is how often callers like Combustor and
    Junction need to) doesn't reintroduce the redundant-mechanism-load cost
    a fresh CanteraFluid would. Composition is fixed at construction (one
    residual evaluation's worth of reacted-product or mixed-stream
    composition); only T/P vary against it, same pattern as CanteraFluid
    itself, just skipping the mechanism (re)load.

    Shared by Combustor (its equilibrium reaction product) and Junction (a
    mass-weighted blend of two or more differently-composed inlet streams)
    -- anywhere a component computes its own fixed composition once per
    residual evaluation and needs a BaseFluid view over it for
    outlet_fluid()/merge_fluids().

    mechanism is carried along (not used internally here) purely so a
    consumer like Junction._blend() can compare it against another fluid's
    mechanism via the shared mass_fractions()/mechanism duck-typed contract
    (see CanteraFluid.mass_fractions()) without needing an isinstance()
    check against this specific class.
    """

    def __init__(self, name: str, gas, mass_fractions: dict[str, float], mechanism: str):
        self.name = name
        self._gas = gas
        self._Y = mass_fractions
        self.mechanism = mechanism

    def mass_fractions(self) -> dict[str, float]:
        # _Y (not self._gas, which is a shared Solution possibly left at
        # some other instance's last-set (T, P) by the time this is read --
        # see this class's own docstring) is always this specific
        # instance's actual, fixed composition.
        return dict(self._Y)

    def _set(self, T: float, P: float) -> None:
        self._gas.TPY = T, P, self._Y

    def enthalpy_pt(self, P: float, T: float) -> float:
        self._set(T, P)
        return float(self._gas.enthalpy_mass)

    def cp(self, P: float, T: float) -> float:
        self._set(T, P)
        return float(self._gas.cp_mass)

    def cv(self, P: float, T: float) -> float:
        self._set(T, P)
        return float(self._gas.cv_mass)

    def temperature_ph(self, P: float, h: float) -> float:
        self._set(_REFERENCE_T, P)
        self._gas.HP = h, P
        return float(self._gas.T)

    def density_ph(self, P: float, h: float) -> float:
        T = self.temperature_ph(P, h)
        self._set(T, P)
        return float(self._gas.density)


class CanteraFluid(BaseFluid):
    """Real-gas mixture fluid backed by a Cantera Solution at a fixed
    composition — temperature-dependent thermo (NASA polynomials) for an
    arbitrary named mixture (humid air, flue gas, natural gas, ...), unlike
    IdealGasFluid/IdealGasMixtureFluid's constant-cp assumption or
    CoolPropFluid's single-species real-fluid tables.

    composition is fixed for the lifetime of this fluid instance — like
    every other fluid here, it's one BaseFluid shared by an entire Network
    (see NetworkState), so this represents a *known, unchanging* mixture's
    real thermo, not a reacting/composition-changing stream. Use Combustor
    (which uses Cantera internally too) to find what a reacting mixture's
    temperature settles at; use this to give the network more accurate
    property lookups for gas already known to be some particular fixed blend
    (e.g. dry vs. humid air, or a representative fixed exhaust composition).

    basis: "mole" (composition is a mole-fraction string/dict, Cantera's
    TPX) or "mass" (mass fractions, TPY). Either way composition can be a
    Cantera composition string ("N2:0.79, O2:0.21") or a {species: fraction}
    dict — anything Cantera's own TPX/TPY setters accept.

    Requires the optional 'cantera' extra: pip install thermowave[cantera]
    """

    def __init__(
        self,
        name: str,
        composition: str | dict[str, float],
        mechanism: str = "gri30.yaml",
        basis: str = "mole",
    ):
        try:
            import cantera as ct
        except ImportError as exc:
            raise ImportError(
                "CanteraFluid requires the 'cantera' extra: pip install thermowave[cantera]"
            ) from exc
        if basis not in ("mole", "mass"):
            raise ValueError(f"CanteraFluid {name!r}: basis must be 'mole' or 'mass', got {basis!r}")

        self.name = name
        self.composition = composition
        self.mechanism = mechanism
        self.basis = basis
        self._gas = ct.Solution(mechanism)
        self._set_composition(_REFERENCE_T, STANDARD_ATMOSPHERE_PA)

    def mass_fractions(self) -> dict[str, float]:
        # Safe to read _gas directly here (unlike _CanteraCompositionFluid's
        # own version of this method): composition is fixed for this
        # class's whole lifetime, and _set_composition() re-asserts that
        # same fixed composition on every property call, so _gas's live
        # mass fractions always match self.composition regardless of
        # whatever (T, P) it was last called at.
        return self._gas.mass_fraction_dict(threshold=0.0)

    def _set_composition(self, T: float, P: float) -> None:
        if self.basis == "mole":
            self._gas.TPX = T, P, self.composition
        else:
            self._gas.TPY = T, P, self.composition

    def enthalpy_pt(self, P: float, T: float) -> float:
        self._set_composition(T, P)
        return float(self._gas.enthalpy_mass)

    def cp(self, P: float, T: float) -> float:
        self._set_composition(T, P)
        return float(self._gas.cp_mass)

    def cv(self, P: float, T: float) -> float:
        self._set_composition(T, P)
        return float(self._gas.cv_mass)

    def temperature_ph(self, P: float, h: float) -> float:
        self._set_composition(_REFERENCE_T, P)
        self._gas.HP = h, P
        return float(self._gas.T)

    def density_ph(self, P: float, h: float) -> float:
        T = self.temperature_ph(P, h)
        self._set_composition(T, P)
        return float(self._gas.density)
