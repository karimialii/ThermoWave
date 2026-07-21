from __future__ import annotations

import math
from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class CheckValve(BaseComponent):
    """One-way flow-restriction valve: lets flow through like Valve in the
    forward direction, and presents a much stiffer resistance to any
    reverse flow — a real boundary condition against backflow, rather than
    just hoping the network's own topology/pressures never ask for it.

    Same K-factor pressure-drop idea as Valve, but with two changes needed
    for a one-way device to make sense at all:

    1. dp = K_eff * rho * v * |v| / 2 (using v*|v| rather than Valve's v**2)
       so the *direction* of the pressure drop tracks the direction of
       flow, not just its magnitude — Valve's own v**2 form is always
       non-negative regardless of which way mdot points, which is fine for
       a symmetric restriction but wrong for a directional one.
    2. K_eff switches between K (forward, v >= 0) and K * reverse_factor
       (reverse, v < 0).

    A real check valve is a hard 0/1 switch (fully open one way, fully shut
    the other); reverse_factor approximates "fully shut" as "resistance
    high enough that any reverse flow the network's pressures could drive
    through it is negligible" rather than an exact zero. That's a
    deliberate trade-off, not a shortcut avoided elsewhere: a hard switch
    would make residuals() non-differentiable exactly at v == 0, which is
    also exactly the operating point most check-valve networks actually
    sit at (fully closed, no flow) — the same reason Junction/Pipe floor a
    division instead of raising ZeroDivisionError at a degenerate but
    physically real state. Flow through the valve scales with
    1/sqrt(K_eff), so reverse_factor=1000.0 (the default) means about 32x
    less reverse flow than the same |dp| would drive forward; raise it
    further if that's still not tight enough for your network's actual
    reverse pressure differential.
    """

    def __init__(
        self, name: str, D: float, K: float, opening: float = 1.0, reverse_factor: float = 1000.0
    ):
        if K <= 0.0:
            raise ValueError(f"CheckValve {name!r}: K must be > 0, got {K}")
        if reverse_factor <= 1.0:
            raise ValueError(
                f"CheckValve {name!r}: reverse_factor must be > 1 (otherwise it isn't "
                f"one-way), got {reverse_factor}"
            )
        self.name = name
        self.D = D
        self.K = K
        self.opening = opening
        self.reverse_factor = reverse_factor
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"
        self._area = math.pi * D**2 / 4

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def residuals(self, state: "NetworkState") -> list[float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        mdot = state.mdot(self._inlet_node)

        rho = state.fluid_at(self._inlet_node).density_ph(P_in, h_in)
        v = mdot / (rho * self._area)
        K = self.K if v >= 0.0 else self.K * self.reverse_factor
        K_eff = K / self.opening**2
        dp_loss = K_eff * (rho * v * abs(v) / 2.0)

        momentum_residual = P_in - P_out - dp_loss
        energy_residual = h_in - h_out
        mass_residual = state.mdot(self._outlet_node) - mdot
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        mdot = state.mdot(self._inlet_node)
        return {
            "power [W]": mdot * (h_out - h_in),
            "PR [-]": P_out / P_in,
            "mdot [kg/s]": mdot,
        }
