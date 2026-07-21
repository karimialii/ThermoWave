from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.core.settings import settings
from thermowave.fluids.two_phase import require_two_phase

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid

_DET_FLOOR = 1.0e-30  # floors the 2x2 mass/energy solve below against a
# singular system (e.g. a fluid model whose density is locally flat in both
# P and h at the current state).


class Drum(BaseComponent):
    """Steam drum: a constant-volume, two-phase reservoir at saturation — a
    boiler's steam/water separator, with genuine liquid-level dynamics.

    Requires a two-phase-capable fluid (CoolPropFluid); checked at residual
    time and at construction.

    Physics: the drum's own (P, h_avg) are two differential states, driven by
    mass and energy conservation on its contents (mass m = fluid.density_ph(P,
    h) * V; specific internal energy u = h - P/rho):
        dm/dt      = sum(mdot_in) - sum(mdot_out)
        d(m*u)/dt  = sum(mdot_in*h_in) - sum(mdot_out*h_out) - heat_loss
    Expanding u = h - P/rho and m/rho = V (constant volume) turns the second
    equation into one relating dP/dt and dh/dt directly; combined with the
    first (via drho/dP, drho/dh, finite-differenced off the fluid model's own
    density_ph() rather than assuming a closed-form equation of state), that's
    two equations in the two unknowns dP/dt, dh/dt, solved as a 2x2 linear
    system in state_derivative(). Each outlet leaves at its own saturation
    enthalpy (h_g for steam, h_f for water), not at the drum's own mixed h.

    Ports (has_riser=True, the natural-circulation-boiler layout):
      feed_in   — feedwater in
      riser_in  — two-phase steam/water return from the evaporator/risers
      steam_out — saturated vapor off the top (h = h_g(P))
      water_out — saturated liquid off the bottom, to the downcomer (h = h_f(P))
    With has_riser=False, riser_in is omitted (a plain feed -> drum -> steam +
    water separator).

    Residuals (4): the two outlet nodes are pinned to the drum pressure and to
    saturated vapor / saturated liquid enthalpy respectively:
      P_steam - P_drum ; h_steam - h_g(P_drum)
      P_water - P_drum ; h_water - h_f(P_drum)
    All mdots are left free, closed externally by whatever the drum connects
    to (a valve, a fixed downstream pressure, ...) — the drum itself doesn't
    dictate how flow splits between steam draw and downcomer circulation.

    Level: intended for normal operation 0 < x < 1 (a partly-full drum). The
    reported level [-] is the liquid volume fraction. The finite-difference
    density partials in state_derivative() can momentarily cross the dome
    boundary if the drum is driven to x≈0 or x≈1 — a known robustness limit.

    A drum's liquid level has no steady-state restoring force: it's a pure
    integrator (net inflow accumulates indefinitely), the same as a real
    drum needs level control to hold a setpoint. Network.solve() is
    therefore singular in the drum's own h — there is no algebraic level to
    solve for. Use Network.solve_transient() instead, which integrates
    (P, h) forward via state_derivative(), exactly the differential-state
    mechanism every dynamic component here shares.

    fluid is used only at construction, to convert level0 (liquid volume
    fraction, default 0.5) into an initial mass quality x0 at P0 and then
    h0 = enthalpy_pq(P0, x0) — differential_parameters() has no access to
    the network's own fluid model. heat_loss [W]: optional constant heat
    loss from the drum's contents (positive = lost).
    """

    def __init__(
        self,
        name: str,
        V: float,
        P0: float,
        fluid: "BaseFluid",
        level0: float = 0.5,
        has_riser: bool = True,
        heat_loss: float | None = None,
    ):
        if V <= 0.0:
            raise ValueError(f"Drum {name!r}: V must be > 0, got {V}")
        if not (0.0 < level0 < 1.0):
            raise ValueError(f"Drum {name!r}: level0 must be in (0, 1), got {level0}")
        require_two_phase(fluid, f"Drum {name!r}")

        self.name = name
        self.V = V
        self.has_riser = has_riser
        self.heat_loss = heat_loss
        self.P0 = settings.pressure_to_si(P0)

        # Seed the differential state: convert liquid volume fraction level0
        # into a mass quality x0 at P0, then to an average enthalpy h0.
        h_f0 = fluid.saturated_liquid_enthalpy(self.P0)
        h_g0 = fluid.saturated_vapor_enthalpy(self.P0)
        rho_f0 = fluid.density_ph(self.P0, h_f0)
        rho_g0 = fluid.density_ph(self.P0, h_g0)
        alpha = level0  # liquid volume fraction
        m_liq = alpha * rho_f0
        m_vap = (1.0 - alpha) * rho_g0
        x0 = m_vap / (m_liq + m_vap)
        self.h0 = h_f0 + x0 * (h_g0 - h_f0)
        # Warm-start guesses for the two outlets (fluid isn't available in
        # guess_outlet_for_pair, so cache the saturated enthalpies here).
        self._h_f0 = h_f0
        self._h_g0 = h_g0

        self._feed_in_node = f"{name}.feed_in"
        self._riser_in_node = f"{name}.riser_in"
        self._steam_out_node = f"{name}.steam_out"
        self._water_out_node = f"{name}.water_out"

    def ports(self) -> dict[str, str]:
        ports = {
            "feed_in": self._feed_in_node,
            "steam_out": self._steam_out_node,
            "water_out": self._water_out_node,
        }
        if self.has_riser:
            ports["riser_in"] = self._riser_in_node
        return ports

    def report_category(self) -> str:
        return "drum"

    def differential_parameters(self) -> dict[str, float]:
        return {"P": self.P0, "h": self.h0}

    def _inlet_nodes(self) -> list[str]:
        nodes = [self._feed_in_node]
        if self.has_riser:
            nodes.append(self._riser_in_node)
        return nodes

    def warm_start_pairs(self) -> list[tuple[str, str]]:
        return [("feed_in", "steam_out"), ("feed_in", "water_out")]

    def guess_outlet_for_pair(
        self, pair: tuple[str, str], P_in: float, h_in: float, mdot: float
    ) -> tuple[float, float]:
        if pair[1] == "steam_out":
            return self.P0, self._h_g0
        if pair[1] == "water_out":
            return self.P0, self._h_f0
        return self.P0, self.h0

    def state_derivative(self, state: "NetworkState") -> dict[str, float]:
        P = state.param(f"{self.name}.P")
        h = state.param(f"{self.name}.h")
        fluid = state.fluid_at(self._feed_in_node)
        rho = fluid.density_ph(P, h)

        eps_P = max(abs(P) * 1.0e-6, 1.0)
        eps_h = max(abs(h) * 1.0e-6, 1.0)
        drho_dP = (fluid.density_ph(P + eps_P, h) - rho) / eps_P
        drho_dh = (fluid.density_ph(P, h + eps_h) - rho) / eps_h

        mdot_steam = state.mdot(self._steam_out_node)
        mdot_water = state.mdot(self._water_out_node)
        _, h_steam = state.node(self._steam_out_node)
        _, h_water = state.node(self._water_out_node)

        mdot_in_total = 0.0
        energy_in = 0.0  # sum of mdot_in*(h_in - h)
        for node in self._inlet_nodes():
            mdot_in = state.mdot(node)
            _, h_in = state.node(node)
            mdot_in_total += mdot_in
            energy_in += mdot_in * (h_in - h)

        mdot_out_total = mdot_steam + mdot_water
        energy_out = mdot_steam * (h_steam - h) + mdot_water * (h_water - h)
        heat_loss = self.heat_loss or 0.0

        # [a11 a12] [dP/dt]   [b1]      mass conservation (dm/dt == mdot_net)
        # [a21 a22] [dh/dt] = [b2]      energy conservation (d(m*u)/dt == ...)
        a11, a12, b1 = self.V * drho_dP, self.V * drho_dh, mdot_in_total - mdot_out_total
        a21, a22, b2 = -self.V, rho * self.V, energy_in - energy_out - heat_loss

        det = a11 * a22 - a12 * a21
        det_safe = det if abs(det) > _DET_FLOOR else _DET_FLOOR
        dPdt = (b1 * a22 - a12 * b2) / det_safe
        dhdt = (a11 * b2 - a21 * b1) / det_safe
        return {"P": dPdt, "h": dhdt}

    def residuals(self, state: "NetworkState") -> list[float]:
        fluid = state.fluid_at(self._feed_in_node)
        require_two_phase(fluid, f"Drum {self.name!r}")
        P_drum = state.param(f"{self.name}.P")

        P_steam, h_steam = state.node(self._steam_out_node)
        P_water, h_water = state.node(self._water_out_node)
        h_g = fluid.saturated_vapor_enthalpy(P_drum)
        h_f = fluid.saturated_liquid_enthalpy(P_drum)

        return [
            P_steam - P_drum,
            h_steam - h_g,
            P_water - P_drum,
            h_water - h_f,
        ]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P = state.param(f"{self.name}.P")
        h = state.param(f"{self.name}.h")
        fluid = state.fluid_at(self._feed_in_node)

        h_f = fluid.saturated_liquid_enthalpy(P)
        h_g = fluid.saturated_vapor_enthalpy(P)
        rho = fluid.density_ph(P, h)
        rho_f = fluid.density_ph(P, h_f)
        # Mass quality by enthalpy (valid in the dome), clamped for reporting.
        x = (h - h_f) / (h_g - h_f) if h_g > h_f else 0.0
        x_rep = min(max(x, 0.0), 1.0)
        v_avg = 1.0 / rho
        v_f = 1.0 / rho_f
        level = (1.0 - x_rep) * v_f / v_avg  # liquid volume fraction

        return {
            "P [Pa]": P,
            "T_sat [K]": fluid.saturation_temperature(P),
            "level [-]": min(max(level, 0.0), 1.0),
            "mdot_steam [kg/s]": state.mdot(self._steam_out_node),
            "mdot_water [kg/s]": state.mdot(self._water_out_node),
            "V [m^3]": self.V,
        }
