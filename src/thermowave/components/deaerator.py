from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.fluids.two_phase import require_two_phase, supports_two_phase

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_MDOT_FLOOR = 1.0e-9  # kg/s, same degenerate-all-zero-inflow guard Junction uses


class Deaerator(BaseComponent):
    """Open feedwater heater: N inlet streams (direct steam extraction,
    condensate drain cascade, ...) mix at one common pressure and leave as a
    single saturated-liquid stream — the open-FWH counterpart to the closed
    FeedwaterHeater, and the physical mixing/equilibrium piece Junction
    alone doesn't provide.

    Junction already does mass-weighted mixing (h_out = mdot-weighted
    average of the inlets) but has no equilibrium constraint at all — a
    Junction standing in for a deaerator today only works if the upstream
    extraction split happens to be sized so the mixed enthalpy lands near
    saturation; nothing checks or corrects it. This component instead pins
    the outlet DIRECTLY to h_f(P) (saturated liquid at the shared outlet
    pressure) — the real physical behavior of an open FWH/deaerator, which
    settles to saturated liquid regardless of exactly what the inlet mix
    averages to (venting off any non-condensables/excess flash steam is
    what makes that true; not modeled explicitly here, just its outcome).

    Same pattern Drum already uses for its own steam_out/water_out pins
    (h_steam - h_g, h_water - h_f), minus Drum's differential-storage
    machinery: this is a purely algebraic, instantaneous mixing+equilibrium
    point, not a dynamic vessel with its own (P, h) state to integrate.

    Common outlet pressure is taken from the first inlet, same simplification
    Junction's own docstring already documents (no attempt to solve for a
    physically self-consistent merge pressure across streams at genuinely
    different pressures).
    """

    def __init__(self, name: str, n_inlets: int):
        if n_inlets < 2:
            raise ValueError(f"Deaerator {name!r}: n_inlets must be >= 2, got {n_inlets}")
        self.name = name
        self.n_inlets = n_inlets
        self._inlet_nodes = [f"{name}.in{i}" for i in range(n_inlets)]
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        ports = {f"in{i}": node for i, node in enumerate(self._inlet_nodes)}
        ports["out"] = self._outlet_node
        return ports

    def warm_start_pairs(self) -> list[tuple[str, str]]:
        # "in"/"out" isn't this component's own port naming (in0/in1/.../out),
        # so BaseComponent's default (which only looks for ports literally
        # named "in"/"out") would silently skip it -- same fix Junction's own
        # warm_start_pairs() override applies, for the same reason.
        if not self._inlet_nodes:
            return []
        return [("in0", "out")]

    def residuals(self, state: "NetworkState") -> list[float]:
        fluid = state.fluid_at(self._inlet_nodes[0])
        require_two_phase(fluid, f"Deaerator {self.name!r}")
        assert supports_two_phase(fluid)

        inlet_mdots = [state.mdot(node) for node in self._inlet_nodes]
        mdot_total = sum(inlet_mdots)
        P_ref, _h_ref = state.node(self._inlet_nodes[0])

        P_out, h_out = state.node(self._outlet_node)
        h_f = fluid.saturated_liquid_enthalpy(P_out)

        momentum_residual = P_out - P_ref
        energy_residual = h_out - h_f
        mass_residual = state.mdot(self._outlet_node) - mdot_total
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_out, h_out = state.node(self._outlet_node)
        fluid = state.fluid_at(self._inlet_nodes[0])
        require_two_phase(fluid, f"Deaerator {self.name!r}")
        assert supports_two_phase(fluid)

        return {
            "P [Pa]": P_out,
            "T_sat [K]": fluid.saturation_temperature(P_out),
            "mdot_out [kg/s]": state.mdot(self._outlet_node),
        }

    # No report_category() override: none of the existing categorized
    # tables' column layouts fit (drum's own has level [-]/V [m^3], neither
    # of which applies to a purely algebraic mixing+equilibrium point with
    # no storage volume) -- report_metrics() above still works standalone
    # (e.g. via SolveResult.state() + component.report_metrics()), just
    # outside the categorized-tables section of print_report().
