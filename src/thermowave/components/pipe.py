from __future__ import annotations

import math
from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_MIN_MDOT = 1.0e-9  # kg/s, floor for heat_loss/mdot below — a Newton iterate
# that clamps a free mdot near zero (Solver.MDOT_MIN) would otherwise blow up
# this division; guess_outlet()'s own mdot argument gets the same floor since
# it's called with the same kind of solver-propagated guess.


class Pipe(BaseComponent):
    """Pipe with Darcy-Weisbach friction pressure drop and optional heat loss.

    Discretized into n_elem sub-elements, each contributing a momentum residual
    (friction pressure drop) and an energy residual (enthalpy balance), all
    computed from the pipe's own (constant, non-branching) inlet mdot. A
    single mass-conservation residual (outlet mdot == inlet mdot) ties the
    outlet's mdot unknown, if free, back to the inlet's.
    """

    def __init__(
        self,
        name: str,
        L: float,
        D: float,
        f: float,
        n_elem: int = 1,
        heat_loss: float | None = None,
    ):
        self.name = name
        self.L = L
        self.D = D
        self.f = f
        self.n_elem = n_elem
        self.heat_loss = heat_loss
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"
        self._area = math.pi * D**2 / 4

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def internal_nodes(self) -> list[str]:
        return [f"{self.name}__mid{i}" for i in range(1, self.n_elem)]

    def _element_nodes(self) -> list[str]:
        return [self._inlet_node, *self.internal_nodes(), self._outlet_node]

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        h_out = h_in if self.heat_loss is None else h_in - self.heat_loss / max(mdot, _MIN_MDOT)
        return P_in, h_out

    def residuals(self, state: "NetworkState") -> list[float]:
        nodes = self._element_nodes()
        elem_L = self.L / self.n_elem
        q_elem = 0.0 if self.heat_loss is None else self.heat_loss / self.n_elem
        mdot = state.mdot(self._inlet_node)
        # One fluid for the whole pipe (its own inlet port's, not each
        # internal element's) -- internal element nodes never get their own
        # node_fluid entry (only port nodes do, see
        # Network._resolve_node_fluid()), and a Pipe doesn't change
        # composition between elements anyway.
        fluid = state.fluid_at(self._inlet_node)

        out: list[float] = []
        for i in range(self.n_elem):
            P_in, h_in = state.node(nodes[i])
            P_out, h_out = state.node(nodes[i + 1])
            rho = fluid.density_ph(P_in, h_in)
            v = mdot / (rho * self._area)
            dp_friction = self.f * (elem_L / self.D) * (rho * v**2 / 2)
            out.append(P_in - P_out - dp_friction)
            out.append(h_in - h_out - q_elem / max(mdot, _MIN_MDOT))
        out.append(state.mdot(self._outlet_node) - mdot)
        return out
