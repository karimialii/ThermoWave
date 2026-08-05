from __future__ import annotations

import math
from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_MIN_OPENING = 1.0e-6  # floor for opening below -- opening=0.0 (fully
# closed, a routine physical state) would otherwise make K_eff = K /
# opening**2 divide by zero; flooring it instead gives a very large but
# finite resistance, letting the solver converge toward (not exactly)
# zero flow, same rationale as Pipe's/Junction's own division floors.


class Valve(BaseComponent):
    """Flow-restriction pressure drop via a loss coefficient K.

    dp = K_eff * (rho * v**2 / 2), where K_eff = K / opening**2 scales the
    resistance up as the valve closes (opening in [0, 1], 1.0 = fully open).
    opening is floored at _MIN_OPENING internally, so opening=0.0 (fully
    closed) gives a very large but finite K_eff instead of raising.
    Throttling is isenthalpic: no work, no heat exchange, so h_in == h_out.
    """

    def __init__(self, name: str, D: float, K: float, opening: float = 1.0):
        if D <= 0.0:
            raise ValueError(f"Valve {name!r}: D must be > 0, got {D}")
        if K < 0.0:
            raise ValueError(f"Valve {name!r}: K must be >= 0, got {K}")
        if not (0.0 <= opening <= 1.0):
            raise ValueError(f"Valve {name!r}: opening must be in [0, 1], got {opening}")
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
        K_eff = self.K / max(self.opening, _MIN_OPENING) ** 2
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
