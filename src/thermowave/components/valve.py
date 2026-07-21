from __future__ import annotations

import math
from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class Valve(BaseComponent):
    """Flow-restriction pressure drop via a loss coefficient K.

    dp = K_eff * (rho * v**2 / 2), where K_eff = K / opening**2 scales the
    resistance up as the valve closes (opening in (0, 1], 1.0 = fully open).
    Throttling is isenthalpic: no work, no heat exchange, so h_in == h_out.
    """

    def __init__(self, name: str, D: float, K: float, opening: float = 1.0):
        self.name = name
        self.D = D
        self.K = K
        self.opening = opening
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
        K_eff = self.K / self.opening**2
        dp_loss = K_eff * (rho * v**2 / 2)

        momentum_residual = P_in - P_out - dp_loss
        energy_residual = h_in - h_out
        mass_residual = state.mdot(self._outlet_node) - mdot
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        return {
            "power [W]": state.mdot(self._inlet_node) * (h_out - h_in),
            "PR [-]": P_out / P_in,
        }
