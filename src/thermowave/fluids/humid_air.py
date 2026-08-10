from thermowave.core.constants import STANDARD_ATMOSPHERE_PA
from thermowave.core.exceptions import FluidRangeError
from thermowave.fluids.base_fluid import BaseFluid

# Dry-air/water-vapor specific gas constants for the cv() approximation below.
# R_DA matches the value ideal_gas_mixture.py's own tests use for air; R_V is
# derived from that same module's _R_UNIVERSAL and H2O's molar mass, rather
# than a second hardcoded constant drifting from it.
R_DA = 287.05  # J/(kg*K), dry air
R_V = 8.314462618 / 0.018015  # J/(kg*K), water vapor (R_universal / M_H2O)


class HumidAirFluid(BaseFluid):
    """Psychrometric (humid-air) fluid model backed by CoolProp's HAPropsSI.

    Requires the optional 'coolprop' extra: pip install thermowave[coolprop]

    W (humidity ratio, kg water / kg dry air) is the fixed, instance-lifetime
    "composition" -- the same role CanteraFluid's composition and
    IdealGasMixtureFluid's mass fractions play, fixed at construction and
    never varying as this fluid instance flows through a network. This is
    the physically correct quantity to fix: a parcel of moist air carries a
    constant W through sensible heating/cooling (only T and P change), while
    relative humidity is *derived* from (T, P, W) and changes continuously
    as T changes -- fixing RH instead would be as wrong as fixing quality
    instead of composition on a CoolPropFluid. Use from_relative_humidity()
    to construct from the more common "50% RH at 20 C" way of specifying an
    intended condition; it derives W once via HAPropsSI and forwards here.

    mdot convention: for a stream carrying this fluid, mdot means DRY-AIR
    mass flow, not total moist-air mass flow -- the ASHRAE/psychrometric
    convention, and the one HAPropsSI's own H/Vha/cp_ha outputs are already
    normalized to (CoolProp exposes 'cp_ha'/'cp' as genuinely different
    per-kg-dry-air vs. per-kg-humid-air keys precisely because psychrometrics
    needs the dry-air basis: dry-air mass is what's actually conserved
    through sensible heating/cooling and through most real HVAC processes
    that don't add/remove water). Building enthalpy_pt/cp/density_ph on the
    _ha-suffixed keys keeps Q = mdot * (h_out - h_in) exact across a
    sensible-only component with mdot held constant. Concretely:
    Source(mdot=5.0, fluid=HumidAirFluid(W=0.01)) means 5 kg/s of dry air
    (5.05 kg/s of actual moist air moving through a duct) -- anyone computing
    duct velocity or comparing against a fan curve needs 1/Vha (moist-air
    volume PER KG DRY AIR) paired with this same dry-air-basis mdot, which is
    exactly what density_ph() below returns; the m^3 cancels out correctly
    because both sides share the dry-air basis.

    No caching: like CoolPropFluid (the closest analog, also a stateless
    free-function wrapper), every property call is a fresh HAPropsSI
    invocation. HAPropsSI has no mutable composition object to amortize the
    way CanteraFluid's Solution does (that class's _composition_set flag
    exists because someone MEASURED a ~180x cost difference for a genuinely
    stateful object re-parsing composition every call -- HAPropsSI has no
    equivalent object to avoid re-touching). If profiling a run with heavy
    HumidAirFluid property-call volume later shows this matters, a
    per-instance memoized (P, T) -> property-vector cache (one HAPropsSI call
    can return several properties at once) would be the natural next step --
    not built here without a demonstrated need, mirroring how CoolPropFluid
    itself ships with no caching either.
    """

    def __init__(
        self,
        name: str,
        W: float,
        P_min: float = 1.0e4,
        P_max: float = 2.0e6,
        T_min: float = 200.0,
        T_max: float = 400.0,
    ):
        try:
            from CoolProp.HumidAirProp import HAPropsSI
        except ImportError as exc:
            raise ImportError(
                "HumidAirFluid requires the 'coolprop' extra: "
                "pip install thermowave[coolprop]"
            ) from exc
        if W < 0.0:
            raise ValueError(f"HumidAirFluid {name!r}: W must be >= 0, got {W}")
        if W > 1.0:
            # W > 1 means more water mass than dry-air mass -- far past any
            # physically reachable saturation at engineering conditions, and
            # well outside HAPropsSI's own correlation range.
            raise ValueError(
                f"HumidAirFluid {name!r}: W={W} is not physically reachable "
                f"(more water mass than dry-air mass); check units (kg/kg, not %)"
            )

        self._ha_props_si = HAPropsSI
        self.name = name
        self.W = W
        self.P_min = P_min
        self.P_max = P_max
        self.T_min = T_min
        self.T_max = T_max

    @classmethod
    def from_relative_humidity(
        cls,
        name: str,
        T: float,
        RH: float,
        P: float = STANDARD_ATMOSPHERE_PA,
        **kwargs,
    ) -> "HumidAirFluid":
        """Construct from a reference (T [K], RH [-, 0..1], P [Pa]) condition
        instead of a direct humidity ratio -- the more common way an intended
        stream condition ("50% RH at 20 C") gets specified. Derives W once
        via HAPropsSI at this reference condition and forwards to __init__;
        further kwargs (P_min/P_max/T_min/T_max) pass through unchanged.
        """
        try:
            from CoolProp.HumidAirProp import HAPropsSI
        except ImportError as exc:
            raise ImportError(
                "HumidAirFluid requires the 'coolprop' extra: "
                "pip install thermowave[coolprop]"
            ) from exc
        try:
            W = HAPropsSI("W", "T", T, "RH", RH, "P", P)
        except ValueError as exc:
            raise FluidRangeError(
                f"HumidAirFluid {name!r}: HAPropsSI W lookup failed for "
                f"T={T}, RH={RH}, P={P}: {exc}"
            ) from exc
        return cls(name, W, **kwargs)

    def _validate_pressure(self, P: float) -> None:
        if P < self.P_min or P > self.P_max:
            raise FluidRangeError(
                f"HumidAirFluid {self.name!r}: pressure {P} outside valid "
                f"range [{self.P_min}, {self.P_max}]"
            )

    def _validate_temperature(self, T: float) -> None:
        if T < self.T_min or T > self.T_max:
            raise FluidRangeError(
                f"HumidAirFluid {self.name!r}: temperature {T} outside valid "
                f"range [{self.T_min}, {self.T_max}]"
            )

    def density_ph(self, P: float, h: float) -> float:
        T = self.temperature_ph(P, h)  # already validates P
        try:
            Vha = self._ha_props_si("Vha", "T", T, "P", P, "W", self.W)
        except ValueError as exc:
            raise FluidRangeError(
                f"HumidAirFluid density_ph failed for P={P}, h={h}, "
                f"fluid={self.name}, W={self.W}: {exc}"
            ) from exc
        return 1.0 / Vha

    def temperature_ph(self, P: float, h: float) -> float:
        self._validate_pressure(P)
        try:
            return self._ha_props_si("T", "P", P, "H", h, "W", self.W)
        except ValueError as exc:
            raise FluidRangeError(
                f"HumidAirFluid temperature_ph failed for P={P}, h={h}, "
                f"fluid={self.name}, W={self.W}: {exc}"
            ) from exc

    def enthalpy_pt(self, P: float, T: float) -> float:
        self._validate_pressure(P)
        self._validate_temperature(T)
        try:
            return self._ha_props_si("H", "T", T, "P", P, "W", self.W)
        except ValueError as exc:
            raise FluidRangeError(
                f"HumidAirFluid enthalpy_pt failed for P={P}, T={T}, "
                f"fluid={self.name}, W={self.W}: {exc}"
            ) from exc

    def cp(self, P: float, T: float) -> float:
        self._validate_pressure(P)
        self._validate_temperature(T)
        try:
            return self._ha_props_si("cp_ha", "T", T, "P", P, "W", self.W)
        except ValueError as exc:
            raise FluidRangeError(
                f"HumidAirFluid cp failed for P={P}, T={T}, "
                f"fluid={self.name}, W={self.W}: {exc}"
            ) from exc

    def cv(self, P: float, T: float) -> float:
        """Ideal-gas approximation: cv = cp - (R_DA + W*R_V).

        HAPropsSI has no native 'cv' output for moist air -- standard, since
        engineering psychrometrics almost never needs it. This treats the
        dry-air + water-vapor mixture as ideal-gas-like for the cp - cv = R
        relation (see module-level R_DA/R_V), which is adequate for
        gamma()-based turbomachinery relations at typical near-atmospheric,
        moderate-temperature HVAC/process conditions, but should not be
        trusted near saturation or at high pressure, where real-gas/
        condensation effects this approximation ignores start to matter --
        use cp()/enthalpy_pt()/density_ph() (all directly HAPropsSI-native)
        for anything precision-sensitive.
        """
        R_eff = R_DA + self.W * R_V
        return self.cp(P, T) - R_eff

    # --- Psychrometric extras (relative humidity / wet-bulb / dew point) ---
    # Additive to HumidAirFluid only, not to the BaseFluid interface (no
    # other fluid model here has a humidity concept at all). Consumers detect
    # these via thermowave.fluids.psychrometrics.supports_humid_air() rather
    # than isinstance -- see that module's docstring.

    def relative_humidity_pt(self, P: float, T: float) -> float:
        """Relative humidity [-, 0..1] at pressure P [Pa], temperature T [K],
        for this fluid's own fixed W."""
        self._validate_pressure(P)
        self._validate_temperature(T)
        try:
            return self._ha_props_si("R", "T", T, "P", P, "W", self.W)
        except ValueError as exc:
            raise FluidRangeError(
                f"HumidAirFluid relative_humidity_pt failed for P={P}, T={T}, "
                f"fluid={self.name}, W={self.W}: {exc}"
            ) from exc

    def wet_bulb_pt(self, P: float, T: float) -> float:
        """Wet-bulb temperature [K] at pressure P [Pa], temperature T [K],
        for this fluid's own fixed W."""
        self._validate_pressure(P)
        self._validate_temperature(T)
        try:
            return self._ha_props_si("B", "T", T, "P", P, "W", self.W)
        except ValueError as exc:
            raise FluidRangeError(
                f"HumidAirFluid wet_bulb_pt failed for P={P}, T={T}, "
                f"fluid={self.name}, W={self.W}: {exc}"
            ) from exc

    def dew_point_pt(self, P: float, T: float) -> float:
        """Dew-point temperature [K] at pressure P [Pa], temperature T [K],
        for this fluid's own fixed W."""
        self._validate_pressure(P)
        self._validate_temperature(T)
        try:
            return self._ha_props_si("D", "T", T, "P", P, "W", self.W)
        except ValueError as exc:
            raise FluidRangeError(
                f"HumidAirFluid dew_point_pt failed for P={P}, T={T}, "
                f"fluid={self.name}, W={self.W}: {exc}"
            ) from exc
