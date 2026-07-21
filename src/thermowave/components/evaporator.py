from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.fluids.two_phase import require_two_phase

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_MIN_MDOT = 1.0e-9  # kg/s, floor for the Q/mdot energy-residual divisions


class Evaporator(BaseComponent):
    """Two-stream evaporator: a boiling working fluid (wf) is heated by an
    explicit single-phase heat source (src) flowing through the same unit —
    a real coupling between two streams, rather than adding heat with no
    accounting for where it came from.

    Four ports: wf_in/wf_out (working fluid, boiling) and src_in/src_out
    (heat source, single-phase). Both streams are tracked through whatever
    fluid model reaches their own inlet node (via NetworkState.fluid_at()) —
    there's no cross-stream property coupling beyond the shared heat duty, so
    this is most physical when both sides are the same real fluid (e.g. hot
    pressurized water boiling lower-pressure water, or a steam-heated
    feedwater evaporator). The wf side must be two-phase-capable
    (CoolPropFluid); checked at residual time.

    Duty comes from the wf outlet spec, never an effectiveness/UA calculation
    — cp is effectively infinite during a constant-pressure phase change, so
    an effectiveness-NTU (C = mdot*cp) approach can't represent boiling. The
    wf outlet is instead pinned to saturated vapor (outlet_quality=1.0) or,
    if superheat > 0, that many K above the saturation temperature:
        Q = mdot_wf * (h_wf_out_target - h_wf_in)
    and the source stream gives up exactly that heat:
        h_src_out = h_src_in - Q / mdot_src

    Pinch is a *diagnostic, not a solved constraint*: report_metrics() exposes
    pinch [K] = T_src_out - T_sat(P_wf_out). Nothing here checks that the
    source can actually supply enough heat at a high enough temperature to
    hit the requested wf outlet spec — the residuals will still solve to a
    self-consistent energy balance even if that requires the source to end up
    colder than the boiling fluid, which is thermodynamically impossible. A
    negative reported pinch is exactly that signal: an infeasible spec, not
    a solver bug.

    Each stream's pressure drop is a fixed ratio (PR_wf/PR_src). 6 residuals:
    wf momentum/energy/mass + src momentum/energy/mass.
    """

    def __init__(
        self,
        name: str,
        PR_wf: float = 1.0,
        PR_src: float = 1.0,
        outlet_quality: float = 1.0,
        superheat: float = 0.0,
    ):
        if not (0.0 < PR_wf <= 1.0):
            raise ValueError(f"Evaporator {name!r}: PR_wf must be in (0, 1], got {PR_wf}")
        if not (0.0 < PR_src <= 1.0):
            raise ValueError(f"Evaporator {name!r}: PR_src must be in (0, 1], got {PR_src}")
        if not (0.0 <= outlet_quality <= 1.0):
            raise ValueError(
                f"Evaporator {name!r}: outlet_quality must be in [0, 1], got {outlet_quality}"
            )
        if superheat < 0.0:
            raise ValueError(f"Evaporator {name!r}: superheat must be >= 0, got {superheat}")
        self.name = name
        self.PR_wf = PR_wf
        self.PR_src = PR_src
        self.outlet_quality = outlet_quality
        self.superheat = superheat
        self._wf_in_node = f"{name}.wf_in"
        self._wf_out_node = f"{name}.wf_out"
        self._src_in_node = f"{name}.src_in"
        self._src_out_node = f"{name}.src_out"

    def ports(self) -> dict[str, str]:
        return {
            "wf_in": self._wf_in_node,
            "wf_out": self._wf_out_node,
            "src_in": self._src_in_node,
            "src_out": self._src_out_node,
        }

    def report_category(self) -> str:
        return "phase_change"

    def warm_start_pairs(self) -> list[tuple[str, str]]:
        return [("wf_in", "wf_out"), ("src_in", "src_out")]

    def guess_outlet_for_pair(
        self, pair: tuple[str, str], P_in: float, h_in: float, mdot: float
    ) -> tuple[float, float]:
        # Crude fixed-swing warm start, only for seeding the Newton solve's
        # initial guess before the real duty is known: the wf side gains a
        # large latent-heat swing boiling up, the source side loses a
        # smaller sensible swing.
        if pair == ("wf_in", "wf_out"):
            return self.PR_wf * P_in, h_in + 2.0e6
        if pair == ("src_in", "src_out"):
            return self.PR_src * P_in, h_in - 3.0e5
        return P_in, h_in

    def _h_wf_out_target(self, fluid, P_wf_out: float) -> float:
        if self.superheat > 0.0:
            return fluid.enthalpy_pt(P_wf_out, fluid.saturation_temperature(P_wf_out) + self.superheat)
        return fluid.enthalpy_pq(P_wf_out, self.outlet_quality)

    def _duty(self, state: "NetworkState") -> float:
        wf_fluid = state.fluid_at(self._wf_in_node)
        _P_wf_in, h_wf_in = state.node(self._wf_in_node)
        P_wf_out, _h_wf_out = state.node(self._wf_out_node)
        mdot_wf = state.mdot(self._wf_in_node)
        h_wf_out_target = self._h_wf_out_target(wf_fluid, P_wf_out)
        return mdot_wf * (h_wf_out_target - h_wf_in)

    def residuals(self, state: "NetworkState") -> list[float]:
        require_two_phase(state.fluid_at(self._wf_in_node), f"Evaporator {self.name!r}")
        wf_fluid = state.fluid_at(self._wf_in_node)

        P_wf_in, _h_wf_in = state.node(self._wf_in_node)
        P_wf_out, h_wf_out = state.node(self._wf_out_node)
        P_src_in, h_src_in = state.node(self._src_in_node)
        P_src_out, h_src_out = state.node(self._src_out_node)
        mdot_wf = state.mdot(self._wf_in_node)
        mdot_src = state.mdot(self._src_in_node)

        h_wf_out_target = self._h_wf_out_target(wf_fluid, P_wf_out)
        Q = self._duty(state)

        wf_momentum = P_wf_out - self.PR_wf * P_wf_in
        wf_energy = h_wf_out - h_wf_out_target
        wf_mass = state.mdot(self._wf_out_node) - mdot_wf

        src_momentum = P_src_out - self.PR_src * P_src_in
        src_energy = h_src_out - (h_src_in - Q / max(mdot_src, _MIN_MDOT))
        src_mass = state.mdot(self._src_out_node) - mdot_src

        return [wf_momentum, wf_energy, wf_mass, src_momentum, src_energy, src_mass]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        wf_fluid = state.fluid_at(self._wf_out_node)
        P_wf_out, h_wf_out = state.node(self._wf_out_node)
        P_src_out, h_src_out = state.node(self._src_out_node)
        src_fluid = state.fluid_at(self._src_out_node)

        T_sat = wf_fluid.saturation_temperature(P_wf_out)
        T_src_out = src_fluid.temperature_ph(P_src_out, h_src_out)
        return {
            "power [W]": self._duty(state),
            "PR [-]": self.PR_wf,
            "x_out [-]": wf_fluid.quality_ph(P_wf_out, h_wf_out),
            "T_sat [K]": T_sat,
            "T_out [K]": wf_fluid.temperature_ph(P_wf_out, h_wf_out),
            "pinch [K]": T_src_out - T_sat,
        }
