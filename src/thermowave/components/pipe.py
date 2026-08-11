from __future__ import annotations

import math
from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.components.heat_transfer import heat_loss_watts
from thermowave.core.constants import STANDARD_GRAVITY_M_S2
from thermowave.core.settings import settings

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_MIN_MDOT = 1.0e-9  # kg/s, floor for heat_loss/mdot below — a Newton iterate
# that clamps a free mdot near zero (Solver.MDOT_MIN) would otherwise blow up
# this division; guess_outlet()'s own mdot argument gets the same floor since
# it's called with the same kind of solver-propagated guess.

_MIN_RE = 1.0e-9  # floor for Reynolds number below, same rationale as _MIN_MDOT
_TRANSITION_RE = 2300.0  # laminar/turbulent switch — a hard switch (small kink
# in the finite-difference Jacobian at this Re) is acceptable here rather than
# a smoothed blend.


class Pipe(BaseComponent):
    """Pipe with Darcy-Weisbach friction pressure drop and optional heat loss.

    Discretized into n_elem sub-elements, each contributing a momentum residual
    (friction pressure drop) and an energy residual (enthalpy balance), all
    computed from the pipe's own (constant, non-branching) inlet mdot. A
    single mass-conservation residual (outlet mdot == inlet mdot) ties the
    outlet's mdot unknown, if free, back to the inlet's.

    The Darcy friction factor is computed from roughness and Reynolds number
    by default (re-evaluated fresh every residual call, since Re depends on
    mdot): f = 64/Re below Re=2300 (laminar), or a roughness-dependent
    turbulent correlation above it — Haaland (explicit) or Colebrook-White
    (implicit fixed-point, seeded from the Haaland value), chosen via
    thermowave.core.settings.friction_correlation ("haaland" by default).
    mu is the fluid's dynamic viscosity [Pa*s], given directly as a constant
    since no BaseFluid model in this codebase exposes viscosity (the same
    trust level Turbine/Compressor's gamma already has).

    f: give it directly to skip the roughness/Re calculation entirely and use
    a fixed friction factor instead — e.g. f=0.0 for a frictionless duct
    (a common way to model a heater/combustor segment via heat_loss alone).
    Leave it None (the default) to use roughness/mu; roughness and mu are
    then both required.

    heat_path: optional Convection/Conduction/Radiation
    (thermowave.components.heat_transfer) for a live-computed loss on top of
    (or instead of) the fixed heat_loss above — e.g. a combustion-chamber
    liner segment radiating/convecting to a ThermalMass wall ring, discretized
    station by station along a chain of single-element Pipes, rather than one
    hand-picked constant. Same accounting as Combustor's own heat_path: both
    heat_loss and heat_path (if given) are summed and split evenly across
    n_elem, same as heat_loss alone always has been. Like Combustor/Turbine/
    Compressor, this component's own endpoint is (self, "in") on whichever
    element the path is meant to represent (n_elem==1 makes that unambiguous;
    with n_elem>1 every element shares the same total loss, split evenly, not
    a per-element pick), so build the path after this Pipe already exists, or
    set it afterwards via pipe.set(heat_path=path).

    z_in/z_out: elevation [m] of the inlet/outlet, both 0.0 by default (no
    effect, exactly today's behavior). Given non-equal, each element also
    gets a hydrostatic/buoyancy term rho*g*(z_out-z_in)/n_elem added to its
    momentum residual alongside dp_friction, using that same element's own
    inlet density -- rising elevation (z_out > z_in) reduces P_out relative
    to P_in the way a real fluid column's weight does, on top of friction;
    falling elevation raises it. This is what lets a closed loop of Pipes
    (a heated riser, a cooled downcomer) circulate under buoyancy alone —
    no pump needed, no new solver machinery: the loop's existing momentum
    residuals already select the equilibrium mdot the same way they always
    do, once density differs between legs enough to matter. See the
    thermosiphon example in docs/examples for a full worked loop.
    """

    def __init__(
        self,
        name: str,
        L: float,
        D: float,
        roughness: float | None = None,
        mu: float | None = None,
        f: float | None = None,
        n_elem: int = 1,
        heat_loss: float | None = None,
        heat_path: BaseComponent | list | tuple | None = None,
        z_in: float = 0.0,
        z_out: float = 0.0,
    ):
        if f is None and (roughness is None or mu is None):
            raise ValueError(
                f"Pipe {name!r}: give f directly, or both roughness and mu "
                "to compute it from Reynolds number"
            )
        self.name = name
        self.L = L
        self.D = D
        self.roughness = roughness
        self.mu = mu
        self.f = f
        self.n_elem = n_elem
        self.heat_loss = heat_loss
        self.heat_path = heat_path
        self.z_in = z_in
        self.z_out = z_out
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"
        self._area = math.pi * D**2 / 4

    def _friction_factor(self, rho: float, v: float) -> float:
        if self.f is not None:
            return self.f
        Re = max(rho * abs(v) * self.D / self.mu, _MIN_RE)
        if Re < _TRANSITION_RE:
            return 64.0 / Re
        rel_roughness = self.roughness / self.D
        x = -1.8 * math.log10((rel_roughness / 3.7) ** 1.11 + 6.9 / Re)
        if settings.resolve_friction_correlation() == "haaland":
            return 1.0 / x**2
        for _ in range(15):
            x = -2.0 * math.log10(rel_roughness / 3.7 + 2.51 * x / Re)
        return 1.0 / x**2

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
        elem_dz = (self.z_out - self.z_in) / self.n_elem
        total_loss = (self.heat_loss or 0.0) + heat_loss_watts(self.heat_path, state)
        q_elem = total_loss / self.n_elem
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
            f = self._friction_factor(rho, v)
            dp_friction = f * (elem_L / self.D) * (rho * v**2 / 2)
            dp_gravity = rho * STANDARD_GRAVITY_M_S2 * elem_dz
            out.append(P_in - P_out - dp_friction - dp_gravity)
            out.append(h_in - h_out - q_elem / max(mdot, _MIN_MDOT))
        out.append(state.mdot(self._outlet_node) - mdot)
        return out
