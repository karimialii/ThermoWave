from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.fluids.two_phase import require_two_phase

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_MIN_MDOT = 1.0e-9  # kg/s, floor for the Q/mdot energy-residual divisions


class Condenser(BaseComponent):
    """Two-stream condenser: a condensing working fluid (wf) is cooled by an
    explicit single-phase coolant (cool) flowing through the same unit.

    Four ports: wf_in/wf_out (working fluid, condensing) and cool_in/cool_out
    (coolant, single-phase). Both streams are tracked through whatever fluid
    model reaches their own inlet node (via NetworkState.fluid_at()) — there's
    no cross-stream property coupling beyond the shared heat duty, so this is
    most physical when both sides are the same real fluid.
    The wf side must be two-phase-capable (CoolPropFluid), checked at
    residual time.

    Duty from the wf outlet spec, never an effectiveness/UA calculation (cp
    is effectively infinite during condensation, so an effectiveness-NTU
    approach can't represent it): the wf outlet is pinned to saturated liquid
    (outlet_quality=0.0) or, if subcool > 0, that many K below saturation:
        Q = mdot_wf * (h_wf_in - h_wf_out_target)   (heat released, >= 0)
        h_cool_out = h_cool_in + Q / mdot_cool

    Pinch is a diagnostic, not a solved constraint: report_metrics() exposes
    pinch [K] = T_sat(P_wf_out) - T_cool_out. Nothing here checks that the
    coolant can actually absorb the requested duty at a low enough
    temperature — a negative reported pinch means the coolant would have to
    end up hotter than the condensing fluid, a thermodynamically infeasible
    spec rather than a solver bug. 6 residuals: wf momentum/energy/mass +
    cool momentum/energy/mass.
    """

    def __init__(
        self,
        name: str,
        PR_wf: float = 1.0,
        PR_cool: float = 1.0,
        outlet_quality: float = 0.0,
        subcool: float = 0.0,
    ):
        if not (0.0 < PR_wf <= 1.0):
            raise ValueError(f"Condenser {name!r}: PR_wf must be in (0, 1], got {PR_wf}")
        if not (0.0 < PR_cool <= 1.0):
            raise ValueError(f"Condenser {name!r}: PR_cool must be in (0, 1], got {PR_cool}")
        if not (0.0 <= outlet_quality <= 1.0):
            raise ValueError(
                f"Condenser {name!r}: outlet_quality must be in [0, 1], got {outlet_quality}"
            )
        if subcool < 0.0:
            raise ValueError(f"Condenser {name!r}: subcool must be >= 0, got {subcool}")
        self.name = name
        self.PR_wf = PR_wf
        self.PR_cool = PR_cool
        self.outlet_quality = outlet_quality
        self.subcool = subcool
        self._wf_in_node = f"{name}.wf_in"
        self._wf_out_node = f"{name}.wf_out"
        self._cool_in_node = f"{name}.cool_in"
        self._cool_out_node = f"{name}.cool_out"

    def ports(self) -> dict[str, str]:
        return {
            "wf_in": self._wf_in_node,
            "wf_out": self._wf_out_node,
            "cool_in": self._cool_in_node,
            "cool_out": self._cool_out_node,
        }

    def report_category(self) -> str:
        return "phase_change"

    def warm_start_pairs(self) -> list[tuple[str, str]]:
        return [("wf_in", "wf_out"), ("cool_in", "cool_out")]

    def guess_outlet_for_pair(
        self, pair: tuple[str, str], P_in: float, h_in: float, mdot: float
    ) -> tuple[float, float]:
        if pair == ("wf_in", "wf_out"):
            return self.PR_wf * P_in, h_in - 2.0e6  # condensing: large enthalpy drop
        if pair == ("cool_in", "cool_out"):
            return self.PR_cool * P_in, h_in + 3.0e5
        return P_in, h_in

    def _h_wf_out_target(self, fluid, P_wf_out: float) -> float:
        if self.subcool > 0.0:
            return fluid.enthalpy_pt(P_wf_out, fluid.saturation_temperature(P_wf_out) - self.subcool)
        return fluid.enthalpy_pq(P_wf_out, self.outlet_quality)

    def _duty(self, state: "NetworkState") -> float:
        wf_fluid = state.fluid_at(self._wf_in_node)
        _P_wf_in, h_wf_in = state.node(self._wf_in_node)
        P_wf_out, _h_wf_out = state.node(self._wf_out_node)
        mdot_wf = state.mdot(self._wf_in_node)
        h_wf_out_target = self._h_wf_out_target(wf_fluid, P_wf_out)
        return mdot_wf * (h_wf_in - h_wf_out_target)

    def residuals(self, state: "NetworkState") -> list[float]:
        require_two_phase(state.fluid_at(self._wf_in_node), f"Condenser {self.name!r}")
        wf_fluid = state.fluid_at(self._wf_in_node)

        P_wf_in, _h_wf_in = state.node(self._wf_in_node)
        P_wf_out, h_wf_out = state.node(self._wf_out_node)
        P_cool_in, h_cool_in = state.node(self._cool_in_node)
        P_cool_out, h_cool_out = state.node(self._cool_out_node)
        mdot_wf = state.mdot(self._wf_in_node)
        mdot_cool = state.mdot(self._cool_in_node)

        h_wf_out_target = self._h_wf_out_target(wf_fluid, P_wf_out)
        Q = self._duty(state)

        wf_momentum = P_wf_out - self.PR_wf * P_wf_in
        wf_energy = h_wf_out - h_wf_out_target
        wf_mass = state.mdot(self._wf_out_node) - mdot_wf

        cool_momentum = P_cool_out - self.PR_cool * P_cool_in
        cool_energy = h_cool_out - (h_cool_in + Q / max(mdot_cool, _MIN_MDOT))
        cool_mass = state.mdot(self._cool_out_node) - mdot_cool

        return [wf_momentum, wf_energy, wf_mass, cool_momentum, cool_energy, cool_mass]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        wf_fluid = state.fluid_at(self._wf_out_node)
        P_wf_out, h_wf_out = state.node(self._wf_out_node)
        P_cool_out, h_cool_out = state.node(self._cool_out_node)
        cool_fluid = state.fluid_at(self._cool_out_node)

        T_sat = wf_fluid.saturation_temperature(P_wf_out)
        T_cool_out = cool_fluid.temperature_ph(P_cool_out, h_cool_out)
        return {
            "power [W]": self._duty(state),  # heat rejected by the wf, positive
            "PR [-]": self.PR_wf,
            "x_out [-]": wf_fluid.quality_ph(P_wf_out, h_wf_out),
            "T_sat [K]": T_sat,
            "T_out [K]": wf_fluid.temperature_ph(P_wf_out, h_wf_out),
            "pinch [K]": T_sat - T_cool_out,
        }
