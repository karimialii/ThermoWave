from abc import ABC, abstractmethod


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
