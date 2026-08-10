from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.fluids.two_phase import require_two_phase, supports_two_phase

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_MIN_MDOT = 1.0e-9  # kg/s, floor for the Q/mdot energy-residual divisions


class FeedwaterHeater(BaseComponent):
    """Closed (shell-and-tube) regenerative feedwater heater: condensing
    extraction steam (hot) heats a single-phase feedwater stream (cold) that
    never mixes with it — the closed-FWH counterpart to Deaerator (the open
    FWH, which DOES mix its streams).

    Structurally identical to Condenser (same 4-port shape, same
    outlet_quality/subcool-target duty calc, same reason it isn't built on
    HeatExchanger's effectiveness-NTU machinery: cp is effectively infinite
    during condensation, so NTU can't represent it) -- a FeedwaterHeater
    genuinely IS a condenser whose "coolant" happens to be feedwater instead
    of ambient cooling water; nothing in the physics distinguishes the two
    cases. What's different is only the framing: TTD (terminal temperature
    difference, T_sat_hot - T_cold_out) is FeedwaterHeater's own name for
    the same quantity Condenser reports as "pinch", and is reported the same
    way here -- a diagnostic computed from Q, not a solved constraint. Given
    a fixed cold-side mdot (matching a real plant's fixed feedwater flow,
    and the design-data-driven input this component was built to replace
    the hand-built SimpleCondenser+SimpleEvaporator+Junction workaround
    for), Q comes out of a genuine SHARED calculation -- both sides computed
    from one hot-side duty, not two independently-fixed numbers that only
    approximately agree. That's the actual value over the 3-component
    workaround: no disclosed hot/cold energy-balance gap, because there's
    only one Q, not two.

    Four ports: hot_in/hot_out (condensing extraction steam) and
    cold_in/cold_out (single-phase feedwater) -- named hot/cold rather than
    wf/cool (unlike Condenser) to match HeatExchanger's own convention and
    exergy.py's default _TWO_STREAM_PORTS lookup, so this plugs into
    exergy_report() for free once added to that module's _TWO_STREAM_HX
    tuple (not done here -- a small, separate follow-up).

    Give at most one of outlet_quality (default 0.0, saturated liquid drain)
    or subcool (> 0, that many K below saturation) -- same spec Condenser
    itself exposes. 6 residuals: hot momentum/energy/mass + cold
    momentum/energy/mass.
    """

    def __init__(
        self,
        name: str,
        PR_hot: float = 1.0,
        PR_cold: float = 1.0,
        outlet_quality: float = 0.0,
        subcool: float = 0.0,
    ):
        if not (0.0 < PR_hot <= 1.0):
            raise ValueError(f"FeedwaterHeater {name!r}: PR_hot must be in (0, 1], got {PR_hot}")
        if not (0.0 < PR_cold <= 1.0):
            raise ValueError(
                f"FeedwaterHeater {name!r}: PR_cold must be in (0, 1], got {PR_cold}"
            )
        if not (0.0 <= outlet_quality <= 1.0):
            raise ValueError(
                f"FeedwaterHeater {name!r}: outlet_quality must be in [0, 1], "
                f"got {outlet_quality}"
            )
        if subcool < 0.0:
            raise ValueError(f"FeedwaterHeater {name!r}: subcool must be >= 0, got {subcool}")
        self.name = name
        self.PR_hot = PR_hot
        self.PR_cold = PR_cold
        self.outlet_quality = outlet_quality
        self.subcool = subcool
        self._hot_in_node = f"{name}.hot_in"
        self._hot_out_node = f"{name}.hot_out"
        self._cold_in_node = f"{name}.cold_in"
        self._cold_out_node = f"{name}.cold_out"

    def ports(self) -> dict[str, str]:
        return {
            "hot_in": self._hot_in_node,
            "hot_out": self._hot_out_node,
            "cold_in": self._cold_in_node,
            "cold_out": self._cold_out_node,
        }

    def report_category(self) -> str:
        return "phase_change"

    def warm_start_pairs(self) -> list[tuple[str, str]]:
        return [("hot_in", "hot_out"), ("cold_in", "cold_out")]

    def guess_outlet_for_pair(
        self, pair: tuple[str, str], P_in: float, h_in: float, mdot: float
    ) -> tuple[float, float]:
        if pair == ("hot_in", "hot_out"):
            return self.PR_hot * P_in, h_in - 2.0e6  # condensing: large enthalpy drop
        if pair == ("cold_in", "cold_out"):
            return self.PR_cold * P_in, h_in + 3.0e5
        return P_in, h_in

    def _h_hot_out_target(self, fluid, P_hot_out: float) -> float:
        if self.subcool > 0.0:
            return fluid.enthalpy_pt(
                P_hot_out, fluid.saturation_temperature(P_hot_out) - self.subcool
            )
        return fluid.enthalpy_pq(P_hot_out, self.outlet_quality)

    def _duty(self, state: "NetworkState") -> float:
        hot_fluid = state.fluid_at(self._hot_in_node)
        _P_hot_in, h_hot_in = state.node(self._hot_in_node)
        P_hot_out, _h_hot_out = state.node(self._hot_out_node)
        mdot_hot = state.mdot(self._hot_in_node)
        h_hot_out_target = self._h_hot_out_target(hot_fluid, P_hot_out)
        return mdot_hot * (h_hot_in - h_hot_out_target)

    def residuals(self, state: "NetworkState") -> list[float]:
        require_two_phase(state.fluid_at(self._hot_in_node), f"FeedwaterHeater {self.name!r}")
        hot_fluid = state.fluid_at(self._hot_in_node)

        P_hot_in, _h_hot_in = state.node(self._hot_in_node)
        P_hot_out, h_hot_out = state.node(self._hot_out_node)
        P_cold_in, h_cold_in = state.node(self._cold_in_node)
        P_cold_out, h_cold_out = state.node(self._cold_out_node)
        mdot_hot = state.mdot(self._hot_in_node)
        mdot_cold = state.mdot(self._cold_in_node)

        h_hot_out_target = self._h_hot_out_target(hot_fluid, P_hot_out)
        Q = self._duty(state)

        hot_momentum = P_hot_out - self.PR_hot * P_hot_in
        hot_energy = h_hot_out - h_hot_out_target
        hot_mass = state.mdot(self._hot_out_node) - mdot_hot

        cold_momentum = P_cold_out - self.PR_cold * P_cold_in
        cold_energy = h_cold_out - (h_cold_in + Q / max(mdot_cold, _MIN_MDOT))
        cold_mass = state.mdot(self._cold_out_node) - mdot_cold

        return [hot_momentum, hot_energy, hot_mass, cold_momentum, cold_energy, cold_mass]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        hot_fluid = state.fluid_at(self._hot_out_node)
        P_hot_out, h_hot_out = state.node(self._hot_out_node)
        P_cold_out, h_cold_out = state.node(self._cold_out_node)
        cold_fluid = state.fluid_at(self._cold_out_node)
        require_two_phase(hot_fluid, f"FeedwaterHeater {self.name!r}")
        assert supports_two_phase(hot_fluid)

        T_sat = hot_fluid.saturation_temperature(P_hot_out)
        T_cold_out = cold_fluid.temperature_ph(P_cold_out, h_cold_out)
        return {
            "power [W]": self._duty(state),  # heat given up by the extraction steam, positive
            "PR [-]": self.PR_hot,
            "x_out [-]": hot_fluid.quality_ph(P_hot_out, h_hot_out),
            "T_sat [K]": T_sat,
            "T_out [K]": hot_fluid.temperature_ph(P_hot_out, h_hot_out),
            "pinch [K]": T_sat - T_cold_out,  # a.k.a. TTD in feedwater-heater terminology
        }
