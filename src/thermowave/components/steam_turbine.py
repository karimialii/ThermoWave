from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.fluids.two_phase import require_entropy

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class SteamTurbine(BaseComponent):
    """Steam turbine: expands a (possibly condensing) vapor via an isentropic
    path scaled by efficiency — the wet-steam-correct counterpart to
    SimpleTurbine, whose ideal-gas T_out = T_in * (1/PR)**((gamma-1)/gamma)
    relation is wrong once the expansion crosses into the two-phase dome.

    Entropy-based: s_in = entropy_ph(P_in, h_in); h_out_isentropic =
    enthalpy_ps(P_out, s_in); the actual work is eta_s * (h_in -
    h_out_isentropic). Requires a fluid exposing entropy_ph / enthalpy_ps
    (CoolPropFluid), checked at residual time.

    Specify exactly one of P_out [Pa] (absolute exhaust pressure) or PR
    (P_in / P_out, > 1). 3 residuals: momentum, energy, mass. report_metrics
    exposes x_out [-] (exhaust quality) — a real steam-turbine concern, since
    excessive wetness erodes the last-stage blades.
    """

    def __init__(
        self,
        name: str,
        P_out: float | None = None,
        PR: float | None = None,
        eta_s: float = 0.85,
    ):
        if (P_out is None) == (PR is None):
            raise ValueError(f"SteamTurbine {name!r}: give exactly one of P_out or PR")
        if P_out is not None and P_out <= 0:
            raise ValueError(f"SteamTurbine {name!r}: P_out must be > 0, got {P_out}")
        if PR is not None and PR <= 1.0:
            raise ValueError(
                f"SteamTurbine {name!r}: PR must be > 1 (a turbine drops pressure), got {PR}"
            )
        if not (0.0 < eta_s <= 1.0):
            raise ValueError(f"SteamTurbine {name!r}: eta_s must be in (0, 1], got {eta_s}")
        self.name = name
        self.P_out = P_out
        self.PR = PR
        self.eta_s = eta_s
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def report_category(self) -> str:
        return "turbomachinery"

    def _P_out_target(self, P_in: float) -> float:
        return self.P_out if self.P_out is not None else P_in / self.PR

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        return self._P_out_target(P_in), h_in

    def residuals(self, state: "NetworkState") -> list[float]:
        fluid = state.fluid_at(self._inlet_node)
        require_entropy(fluid, f"SteamTurbine {self.name!r}")
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)

        s_in = fluid.entropy_ph(P_in, h_in)
        h_out_isentropic = fluid.enthalpy_ps(P_out, s_in)
        dh_actual = self.eta_s * (h_in - h_out_isentropic)

        momentum_residual = P_out - self._P_out_target(P_in)
        energy_residual = h_out - (h_in - dh_actual)
        mass_residual = state.mdot(self._outlet_node) - state.mdot(self._inlet_node)
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        fluid = state.fluid_at(self._outlet_node)
        return {
            "power [W]": state.mdot(self._inlet_node) * (h_in - h_out),
            "eta_s [-]": self.eta_s,
            "PR [-]": P_in / P_out,
            "x_out [-]": fluid.quality_ph(P_out, h_out),
        }
