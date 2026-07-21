from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.fluids.two_phase import require_two_phase

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_MIN_MDOT = 1.0e-9  # kg/s, floor for the duty/mdot division below
_H_SWING_GUESS = 2.0e6  # J/kg, crude warm-start enthalpy rise (boiling is a big swing)


class SimpleEvaporator(BaseComponent):
    """Single-stream evaporator / boiler / superheater: adds heat to boil the
    working fluid up to a specified outlet condition, without modeling the
    heat source explicitly (the phase-change counterpart to SimpleCombustor,
    which adds Q without modeling a second stream).

    Requires a two-phase-capable fluid (CoolPropFluid) — the ideal-gas and
    Cantera models have no saturation dome. Checked at residual time via
    thermowave.fluids.two_phase.require_two_phase().

    Duty is NOT computed from an effectiveness/UA calc (cp -> infinity in the
    two-phase dome breaks that framework) — it comes from the outlet state:

    - Outlet-spec mode (default, duty=None): the outlet is pinned to a target
      set by outlet_quality (0..1; 1.0 = saturated vapor) and, if superheat >
      0, that many kelvin above the saturation temperature (a superheater):
        superheat == 0:  h_out = enthalpy_pq(P_out, outlet_quality)
        superheat  > 0:  h_out = enthalpy_pt(P_out, T_sat(P_out) + superheat)
      The heat added Q = mdot * (h_out - h_in) is then a *reported* result.

    - Duty mode (duty given, [W]): the heat added is fixed and the outlet
      enthalpy follows, h_out = h_in + duty / mdot; the resulting outlet
      quality / superheat is reported.

    P_out = PR * P_in (PR default 1.0 — boilers are ~isobaric; set PR < 1 for
    a pressure drop). Mass is conserved (no fuel added), so 3 residuals:
    momentum, energy, mass.
    """

    def __init__(
        self,
        name: str,
        PR: float = 1.0,
        outlet_quality: float = 1.0,
        superheat: float = 0.0,
        duty: float | None = None,
    ):
        if not (0.0 < PR <= 1.0):
            raise ValueError(f"SimpleEvaporator {name!r}: PR must be in (0, 1], got {PR}")
        if not (0.0 <= outlet_quality <= 1.0):
            raise ValueError(
                f"SimpleEvaporator {name!r}: outlet_quality must be in [0, 1], got {outlet_quality}"
            )
        if superheat < 0.0:
            raise ValueError(
                f"SimpleEvaporator {name!r}: superheat must be >= 0, got {superheat}"
            )
        self.name = name
        self.PR = PR
        self.outlet_quality = outlet_quality
        self.superheat = superheat
        self.duty = duty
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def report_category(self) -> str:
        return "phase_change"

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        if self.duty is not None:
            return self.PR * P_in, h_in + self.duty / max(mdot, _MIN_MDOT)
        return self.PR * P_in, h_in + _H_SWING_GUESS

    def _h_out_target(self, state: "NetworkState", P_out: float, h_in: float, mdot: float) -> float:
        fluid = state.fluid_at(self._inlet_node)
        if self.duty is not None:
            return h_in + self.duty / max(mdot, _MIN_MDOT)
        if self.superheat > 0.0:
            T_target = fluid.saturation_temperature(P_out) + self.superheat
            return fluid.enthalpy_pt(P_out, T_target)
        return fluid.enthalpy_pq(P_out, self.outlet_quality)

    def residuals(self, state: "NetworkState") -> list[float]:
        require_two_phase(state.fluid_at(self._inlet_node), f"SimpleEvaporator {self.name!r}")
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        mdot_in = state.mdot(self._inlet_node)

        h_out_target = self._h_out_target(state, P_out, h_in, mdot_in)

        momentum_residual = P_out - self.PR * P_in
        energy_residual = h_out - h_out_target
        mass_residual = state.mdot(self._outlet_node) - mdot_in
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        mdot = state.mdot(self._inlet_node)
        fluid = state.fluid_at(self._outlet_node)

        T_out = fluid.temperature_ph(P_out, h_out)
        T_sat = fluid.saturation_temperature(P_out)
        x_out = fluid.quality_ph(P_out, h_out)
        return {
            "power [W]": mdot * (h_out - h_in),
            "PR [-]": P_out / P_in,
            "x_out [-]": x_out,
            "T_sat [K]": T_sat,
            "T_out [K]": T_out,
            "dT_sat [K]": T_out - T_sat,  # superheat (positive when above saturation)
        }
