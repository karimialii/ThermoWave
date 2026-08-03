import math
from abc import ABC, abstractmethod

# Reference state for ConstantCpFluid's entropy_ph/enthalpy_ps below: entropy
# needs its own finite datum (s=0 at T_ref, P_ref) separate from enthalpy's
# own h=0-at-T=0 datum, since ideal-gas entropy diverges as T->0 or P->0.
# Only entropy DIFFERENCES ever matter physically, so the exact values here
# are arbitrary but fixed.
_S_REF_T = 298.15  # K
_S_REF_P = 101325.0  # Pa
_S_REF = 0.0  # J/(kg*K)


class BaseFluid(ABC):
    """Interface every fluid property model implements.

    All arguments and return values are strict SI units: Pa, K, J/kg, J/(kg*K).
    """

    name: str

    @abstractmethod
    def density_ph(self, P: float, h: float) -> float:
        """Density [kg/m^3] given pressure [Pa] and specific enthalpy [J/kg]."""

    @abstractmethod
    def temperature_ph(self, P: float, h: float) -> float:
        """Temperature [K] given pressure [Pa] and specific enthalpy [J/kg]."""

    @abstractmethod
    def enthalpy_pt(self, P: float, T: float) -> float:
        """Specific enthalpy [J/kg] given pressure [Pa] and temperature [K]."""

    @abstractmethod
    def cp(self, P: float, T: float) -> float:
        """Specific heat at constant pressure [J/(kg*K)]."""

    @abstractmethod
    def cv(self, P: float, T: float) -> float:
        """Specific heat at constant volume [J/(kg*K)]."""

    def gamma(self, P: float, T: float) -> float:
        """Ratio of specific heats cp/cv [-] — used by map-based/analytic
        turbomachinery components (Compressor, Turbine, SimpleCompressor,
        SimpleTurbine) for their isentropic-relation math when they aren't
        given a fixed gamma directly. One implementation here (cp/cv) rather
        than duplicated per fluid model; only cv() varies by model.
        """
        return self.cp(P, T) / self.cv(P, T)


class ConstantCpFluid(BaseFluid):
    """Calorically-perfect gas: constant cp, h = cp * T referenced to h=0 at T=0.

    Subclasses must set self.R [J/(kg*K)] and self._cp [J/(kg*K)] in __init__.
    """

    R: float
    _cp: float

    def cp(self, P: float, T: float) -> float:
        return self._cp

    def cv(self, P: float, T: float) -> float:
        # Ideal-gas relation: cp - cv = R.
        return self._cp - self.R

    def enthalpy_pt(self, P: float, T: float) -> float:
        return self._cp * T

    def temperature_ph(self, P: float, h: float) -> float:
        return h / self._cp

    def density_ph(self, P: float, h: float) -> float:
        T = self.temperature_ph(P, h)
        return P / (self.R * T)

    def entropy_ph(self, P: float, h: float) -> float:
        """Specific entropy [J/(kg*K)], closed-form ideal-gas relation
        referenced to s=0 at (_S_REF_T, _S_REF_P) -- see this module's own
        comment for why that datum is separate from enthalpy_pt's h=0-at-T=0
        one. Gives every ConstantCpFluid subclass (IdealGasFluid,
        IdealGasMixtureFluid) entropy/enthalpy_ps for free, which in turn
        lets Pump/SteamTurbine's entropy-based isentropic path (and any
        exergy calculation) work with a constant-cp ideal gas, not only
        CoolPropFluid -- see thermowave.fluids.two_phase.require_entropy(),
        which duck-types on these two methods existing at all.
        """
        T = self.temperature_ph(P, h)
        return self._cp * math.log(T / _S_REF_T) - self.R * math.log(P / _S_REF_P) + _S_REF

    def enthalpy_ps(self, P: float, s: float) -> float:
        """Inverse of entropy_ph at fixed P: solve for T, then reuse
        enthalpy_pt (not a re-derived cp*T) so a future subclass override of
        enthalpy_pt stays consistent with this."""
        T = _S_REF_T * math.exp((s - _S_REF + self.R * math.log(P / _S_REF_P)) / self._cp)
        return self.enthalpy_pt(P, T)
