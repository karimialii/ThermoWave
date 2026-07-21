from thermowave.core.exceptions import FluidRangeError
from thermowave.fluids.base_fluid import BaseFluid


class CoolPropFluid(BaseFluid):
    """Real-fluid model backed by CoolProp's PropsSI.

    Requires the optional 'coolprop' extra: pip install thermowave[coolprop]
    """

    def __init__(self, name: str, P_min: float = 1.0e3, P_max: float = 1.0e8):
        try:
            from CoolProp.CoolProp import PropsSI
        except ImportError as exc:
            raise ImportError(
                "CoolPropFluid requires the 'coolprop' extra: "
                "pip install thermowave[coolprop]"
            ) from exc

        self._props_si = PropsSI
        self.name = name
        self.P_min = P_min
        self.P_max = P_max

    def _clamp_pressure(self, P: float) -> float:
        return min(max(P, self.P_min), self.P_max)

    def density_ph(self, P: float, h: float) -> float:
        P = self._clamp_pressure(P)
        try:
            return self._props_si("D", "P", P, "H", h, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp density_ph failed for P={P}, h={h}, fluid={self.name}: {exc}"
            ) from exc

    def temperature_ph(self, P: float, h: float) -> float:
        P = self._clamp_pressure(P)
        try:
            return self._props_si("T", "P", P, "H", h, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp temperature_ph failed for P={P}, h={h}, fluid={self.name}: {exc}"
            ) from exc

    def enthalpy_pt(self, P: float, T: float) -> float:
        P = self._clamp_pressure(P)
        try:
            return self._props_si("H", "P", P, "T", T, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp enthalpy_pt failed for P={P}, T={T}, fluid={self.name}: {exc}"
            ) from exc

    def cp(self, P: float, T: float) -> float:
        P = self._clamp_pressure(P)
        try:
            return self._props_si("C", "P", P, "T", T, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp cp failed for P={P}, T={T}, fluid={self.name}: {exc}"
            ) from exc

    def cv(self, P: float, T: float) -> float:
        P = self._clamp_pressure(P)
        try:
            return self._props_si("O", "P", P, "T", T, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp cv failed for P={P}, T={T}, fluid={self.name}: {exc}"
            ) from exc

    # --- Saturation / quality (two-phase) -----------------------------------
    # Additive to CoolPropFluid only, not to the BaseFluid interface (the
    # ideal-gas/Cantera models can't implement them). Consumers detect these
    # via thermowave.fluids.two_phase.supports_two_phase() rather than
    # isinstance -- see that module's docstring.

    def saturation_temperature(self, P: float) -> float:
        """Saturation (boiling) temperature [K] at pressure P [Pa]."""
        P = self._clamp_pressure(P)
        try:
            return self._props_si("T", "P", P, "Q", 0, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp saturation_temperature failed for P={P}, fluid={self.name}: {exc}"
            ) from exc

    def saturation_pressure(self, T: float) -> float:
        """Saturation (vapor) pressure [Pa] at temperature T [K]."""
        try:
            return self._props_si("P", "T", T, "Q", 0, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp saturation_pressure failed for T={T}, fluid={self.name}: {exc}"
            ) from exc

    def saturated_liquid_enthalpy(self, P: float) -> float:
        """Saturated-liquid specific enthalpy h_f [J/kg] at pressure P [Pa]."""
        P = self._clamp_pressure(P)
        try:
            return self._props_si("H", "P", P, "Q", 0, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp saturated_liquid_enthalpy failed for P={P}, fluid={self.name}: {exc}"
            ) from exc

    def saturated_vapor_enthalpy(self, P: float) -> float:
        """Saturated-vapor specific enthalpy h_g [J/kg] at pressure P [Pa]."""
        P = self._clamp_pressure(P)
        try:
            return self._props_si("H", "P", P, "Q", 1, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp saturated_vapor_enthalpy failed for P={P}, fluid={self.name}: {exc}"
            ) from exc

    def enthalpy_pq(self, P: float, x: float) -> float:
        """Specific enthalpy [J/kg] at pressure P [Pa] and vapor quality x [-]
        (0 = saturated liquid, 1 = saturated vapor)."""
        P = self._clamp_pressure(P)
        try:
            return self._props_si("H", "P", P, "Q", x, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp enthalpy_pq failed for P={P}, x={x}, fluid={self.name}: {exc}"
            ) from exc

    def quality_ph(self, P: float, h: float) -> float:
        """Vapor quality x [-] at pressure P [Pa] and enthalpy h [J/kg].

        Inside the two-phase dome this is 0..1. CoolProp returns -1.0 for a
        single-phase state -- BOTH subcooled liquid and superheated vapor
        (it does not return >1 for superheat) -- so callers that need to
        distinguish subcooled from superheated should compare temperature
        against saturation_temperature() rather than relying on the sign
        here.
        """
        P = self._clamp_pressure(P)
        try:
            return self._props_si("Q", "P", P, "H", h, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp quality_ph failed for P={P}, h={h}, fluid={self.name}: {exc}"
            ) from exc

    # --- Entropy (isentropic-path components: Pump, SteamTurbine) ------------

    def entropy_ph(self, P: float, h: float) -> float:
        """Specific entropy [J/(kg*K)] at pressure P [Pa], enthalpy h [J/kg]."""
        P = self._clamp_pressure(P)
        try:
            return self._props_si("S", "P", P, "H", h, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp entropy_ph failed for P={P}, h={h}, fluid={self.name}: {exc}"
            ) from exc

    def enthalpy_ps(self, P: float, s: float) -> float:
        """Specific enthalpy [J/kg] at pressure P [Pa], entropy s [J/(kg*K)]."""
        P = self._clamp_pressure(P)
        try:
            return self._props_si("H", "P", P, "S", s, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp enthalpy_ps failed for P={P}, s={s}, fluid={self.name}: {exc}"
            ) from exc
