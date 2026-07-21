from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.core.settings import settings

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid

_DET_FLOOR = 1.0e-30  # guards the 2x2 solve below against a genuinely
# singular mass/energy system (e.g. a fluid model whose density is locally
# flat in both P and h at the current state) the same way Junction/Pipe
# floor a division against a degenerate Newton iterate.


class Tank(BaseComponent):
    """Constant-volume plenum with real mass/energy storage — the genuinely
    transient counterpart to Junction's zero-volume, quasi-steady mixing.

    Every other component in this network is quasi-steady: its outlet state
    is an algebraic function of its inlet state at the current instant, with
    no memory of its own. A Tank actually accumulates — its contents' own
    (P, h) lag behind whatever's flowing in and out, exactly like a real
    plenum, receiver, or accumulator volume. That's what lets a network
    show genuine pressure/temperature transients (filling, blow-down, surge
    behavior) instead of only ever being a sequence of independent steady
    states — see BaseComponent.differential_parameters()'s docstring for
    the general mechanism (the same one a dynamic Shaft's rotor speed uses).

    Physics: the tank's own (P, h) are its two differential states, driven
    by mass and energy conservation on its contents (mass m = fluid.
    density_ph(P, h) * V; specific internal energy u = h - P/rho):
        dm/dt          = mdot_in - mdot_out
        d(m*u)/dt      = mdot_in*h_in - mdot_out*h_out - heat_loss
    Expanding u = h - P/rho and m/rho = V (constant volume) turns the second
    equation into one relating dP/dt and dh/dt directly; combined with the
    first (via drho/dP, drho/dh — finite-differenced off the fluid model's
    own density_ph(), the same finite-difference philosophy Solver's own
    Jacobian uses, rather than assuming a closed-form ideal-gas equation of
    state), that's two equations in the two unknowns dP/dt, dh/dt, solved as
    a plain 2x2 linear system in state_derivative() below. Works with any
    BaseFluid (ideal gas, CoolProp, Cantera), not just calorically-perfect
    gases.

    Single inlet, single outlet, both left entirely free: the outlet always
    carries the tank's own current (well-mixed) state (residuals() below),
    but — unlike Junction, which fixes outlet mdot as a split fraction of
    inflow — neither mdot is constrained by the tank itself. Real inflow and
    outflow generally differ (that's the entire point of a control volume
    with accumulation); each is instead closed by whatever's upstream/
    downstream (a Source, a Valve, a fixed Sink pressure, ...), the same way
    a Source(mdot=None) relies on something else to pin its own flow.

    Network.solve_transient() integrates (P, h) forward in time (backward-
    Euler, like every other differential state). Network.solve() (steady
    state) instead solves for the (P, h) at which both derivatives are
    exactly zero — a genuine steady-state plenum condition (mdot_in ==
    mdot_out, and h at whatever value keeps d(m*u)/dt balanced), not a
    placeholder.

    fluid is only used once, at construction, to convert P0/T0 into the
    differential state's initial (P, h) value — differential_parameters()
    has no access to the network's fluid model itself (see its docstring),
    unlike guess_free_parameters()), so it's supplied directly here instead;
    pass the same fluid instance given to Network(fluid=...).

    heat_loss [W]: optional constant heat loss from the tank's contents,
    same sign convention as Pipe's own heat_loss (positive = heat lost).
    """

    def __init__(
        self,
        name: str,
        V: float,
        P0: float,
        T0: float,
        fluid: "BaseFluid",
        heat_loss: float | None = None,
    ):
        if V <= 0.0:
            raise ValueError(f"Tank {name!r}: V must be > 0, got {V}")

        self.name = name
        self.V = V
        self.heat_loss = heat_loss
        self.P0 = settings.pressure_to_si(P0)
        T0_si = settings.temperature_to_si(T0)
        self.h0 = fluid.enthalpy_pt(self.P0, T0_si)
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def report_category(self) -> str:
        return "tank"

    def differential_parameters(self) -> dict[str, float]:
        return {"P": self.P0, "h": self.h0}

    def state_derivative(self, state: "NetworkState") -> dict[str, float]:
        P = state.param(f"{self.name}.P")
        h = state.param(f"{self.name}.h")
        fluid = state.fluid_at(self._inlet_node)
        rho = fluid.density_ph(P, h)

        eps_P = max(abs(P) * 1.0e-6, 1.0)
        eps_h = max(abs(h) * 1.0e-6, 1.0)
        drho_dP = (fluid.density_ph(P + eps_P, h) - rho) / eps_P
        drho_dh = (fluid.density_ph(P, h + eps_h) - rho) / eps_h

        mdot_in = state.mdot(self._inlet_node)
        _, h_in = state.node(self._inlet_node)
        mdot_out = state.mdot(self._outlet_node)
        heat_loss = self.heat_loss or 0.0

        # [a11 a12] [dP/dt]   [b1]      mass conservation (dm/dt == mdot_net)
        # [a21 a22] [dh/dt] = [b2]      energy conservation (d(m*u)/dt == ...)
        a11, a12, b1 = self.V * drho_dP, self.V * drho_dh, mdot_in - mdot_out
        a21, a22, b2 = -self.V, rho * self.V, mdot_in * (h_in - h) - heat_loss

        det = a11 * a22 - a12 * a21
        det_safe = det if abs(det) > _DET_FLOOR else _DET_FLOOR
        dPdt = (b1 * a22 - a12 * b2) / det_safe
        dhdt = (a11 * b2 - a21 * b1) / det_safe
        return {"P": dPdt, "h": dhdt}

    def residuals(self, state: "NetworkState") -> list[float]:
        P_tank = state.param(f"{self.name}.P")
        h_tank = state.param(f"{self.name}.h")
        P_out, h_out = state.node(self._outlet_node)
        return [P_out - P_tank, h_out - h_tank]

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        return self.P0, self.h0

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P = state.param(f"{self.name}.P")
        h = state.param(f"{self.name}.h")
        return {
            "P [Pa]": P,
            "T [K]": state.fluid_at(self._inlet_node).temperature_ph(P, h),
            "mdot_in [kg/s]": state.mdot(self._inlet_node),
            "mdot_out [kg/s]": state.mdot(self._outlet_node),
            "V [m^3]": self.V,
        }
