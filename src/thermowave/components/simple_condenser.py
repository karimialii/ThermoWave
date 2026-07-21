from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.fluids.two_phase import require_two_phase

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_MIN_MDOT = 1.0e-9  # kg/s, floor for the duty/mdot division below
_H_SWING_GUESS = 2.0e6  # J/kg, crude warm-start enthalpy drop (condensing is a big swing)


class SimpleCondenser(BaseComponent):
    """Single-stream condenser: rejects heat to condense the working fluid
    down to a specified outlet condition, without modeling the coolant that
    absorbs that heat.

    Requires a two-phase-capable fluid (CoolPropFluid); checked at residual
    time via thermowave.fluids.two_phase.require_two_phase().

    Duty is NOT computed from an effectiveness/UA calculation — cp is
    effectively infinite during a constant-pressure phase change, so an
    effectiveness-NTU (C = mdot*cp) framework can't represent condensation.
    It comes from the outlet state instead:

    - Outlet-spec mode (default, duty=None): outlet pinned to a target set by
      outlet_quality (0..1; 0.0 = saturated liquid) and, if subcool > 0, that
      many kelvin BELOW the saturation temperature (a subcooler):
        subcool == 0:  h_out = enthalpy_pq(P_out, outlet_quality)
        subcool  > 0:  h_out = enthalpy_pt(P_out, T_sat(P_out) - subcool)
      The heat rejected Q = mdot * (h_out - h_in) (negative) is reported.

    - Duty mode (duty given, [W], should be negative for heat rejection):
      h_out = h_in + duty / mdot; resulting outlet quality/subcool reported.

    P_out = PR * P_in (PR default 1.0). Mass conserved; 3 residuals.
    """

    def __init__(
        self,
        name: str,
        PR: float = 1.0,
        outlet_quality: float = 0.0,
        subcool: float = 0.0,
        duty: float | None = None,
    ):
        if not (0.0 < PR <= 1.0):
            raise ValueError(f"SimpleCondenser {name!r}: PR must be in (0, 1], got {PR}")
        if not (0.0 <= outlet_quality <= 1.0):
            raise ValueError(
                f"SimpleCondenser {name!r}: outlet_quality must be in [0, 1], got {outlet_quality}"
            )
        if subcool < 0.0:
            raise ValueError(f"SimpleCondenser {name!r}: subcool must be >= 0, got {subcool}")
        self.name = name
        self.PR = PR
        self.outlet_quality = outlet_quality
        self.subcool = subcool
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
        return self.PR * P_in, h_in - _H_SWING_GUESS

    def _h_out_target(self, state: "NetworkState", P_out: float, h_in: float, mdot: float) -> float:
        fluid = state.fluid_at(self._inlet_node)
        if self.duty is not None:
            return h_in + self.duty / max(mdot, _MIN_MDOT)
        if self.subcool > 0.0:
            T_target = fluid.saturation_temperature(P_out) - self.subcool
            return fluid.enthalpy_pt(P_out, T_target)
        return fluid.enthalpy_pq(P_out, self.outlet_quality)

    def residuals(self, state: "NetworkState") -> list[float]:
        require_two_phase(state.fluid_at(self._inlet_node), f"SimpleCondenser {self.name!r}")
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
            "power [W]": mdot * (h_out - h_in),  # negative: heat rejected
            "PR [-]": P_out / P_in,
            "x_out [-]": x_out,
            "T_sat [K]": T_sat,
            "T_out [K]": T_out,
            "dT_sat [K]": T_sat - T_out,  # subcool (positive when below saturation)
        }
