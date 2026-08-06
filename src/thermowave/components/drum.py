from __future__ import annotations

from typing import TYPE_CHECKING, Callable

from thermowave.components.base_component import BaseComponent
from thermowave.core.exceptions import FluidRangeError
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
    reported level [-] is the liquid volume fraction. state_derivative()'s
    finite-difference density partials are stepped one-sided so they never
    cross the dome boundary even at x≈0 or x≈1 — see that method for why a
    straddling step is not merely inaccurate but sign-flipping.

    Usable in an ordinary Network.solve() (dt=None), not just
    solve_transient(): state_derivative()==0 there reduces exactly to the
    drum's own mass and energy balance (the same two equations a purely
    algebraic, steady-state drum model -- e.g. TESPy's own Drum -- computes
    directly), not something that needs the dynamic-state machinery running.

    What it does NOT do on its own, in EITHER solve mode, is anchor its own
    (P, h): a bare Drum contributes exactly 4 port-pinning residuals (P/h at
    steam_out and water_out relative to whatever P_drum happens to be) plus
    those 2 balance equations -- 6 equations for what works out to (once its
    neighbors' own residuals are accounted for) one MORE free unknown than
    that closes: 3 quantities with no natural closing equation anywhere in a
    typical surrounding topology (this component's own pressure, and, with
    has_riser=True, the riser recirculation mass flow) once oil/HTF mdot on
    a coupled two-stream evaporator is ALSO left free -- traced by hand
    while building the SEGS boiler benchmark (see segs_exergy_benchmark.py's
    BOILER SUBSYSTEM), where neither Setpoint nor Controller can reach a raw
    node mdot or a differential parameter to close the gap generically.
    P_target and water_out_mdot below are the narrow, physically-grounded
    fix: a real boiler drum IS essentially isobaric with its feed (no
    Setpoint/Controller detour needed to say so directly), and the riser
    recirculation ratio is a genuine design input in practice (the same
    role TESPy's own eva.set_attr(ttd_l=5) plays) -- not something a
    conservation-only model derives on its own regardless of the tool.

    P_target [Pa]: when given, adds a 5th residual pinning the drum's own
    differential-state pressure directly (P_drum - P_target), instead of
    leaving it a genuinely free unknown with nothing else in a typical
    topology to anchor it. Leave at the default None for a network that
    already anchors this drum's pressure some other way (or for
    solve_transient(), where dt drives P forward from P0 instead).

    water_out_mdot [kg/s]: when given (has_riser=True only), adds a 6th
    residual pinning water_out's own mass flow directly (the riser/downcomer
    supply rate) -- the recirculation-ratio equivalent of Source's own
    fixed mdot, for exactly the same reason Source needs one: nothing else
    conserves mass INTO a specific numeric value on its own. Leave at the
    default None when something else in the topology already closes it
    (e.g. a fixed-mdot boundary on that branch) or for solve_transient().

    level_target [-]: when given, adds a residual pinning the drum's liquid
    volume fraction (the same quantity report_metrics() reports as
    "level [-]") — the steady-state closure for the drum's OWN differential
    enthalpy state, `h`. This one is not merely convenient like P_target and
    water_out_mdot are; without it a steady solve (dt=None) is EXACTLY
    singular, for a reason worth spelling out because the symptom looks
    nothing like the cause. In steady state the solver closes each
    differential parameter with state_derivative()==0, i.e. this drum's own
    mass and energy balance:
        b1 = sum(mdot_in) - sum(mdot_out)
        b2 = sum(mdot_in*(h_in - h)) - sum(mdot_out*(h_out - h)) - heat_loss
    b1 has no `h` in it at all, and b2's own `h` terms collapse to -h*b1 —
    so d(b2)/dh = -b1, which is identically zero at any solution. `h`
    appears in no other residual either (residuals() pins the two outlets to
    h_g/h_f of P_drum, not of h). The whole Jacobian column for this
    component's `h` therefore vanishes at the solution: a steady drum's
    water inventory/level is genuinely undetermined by conservation alone,
    exactly as a real boiler drum's is (which is why real ones have level
    control). state_derivative()'s 2x2 solve divides through by det, which
    hides this — instead of a clean singular-matrix error, Newton gets a
    tiny-but-nonzero column, lets `h` drift, walks it across h_f(P) out of
    the two-phase dome (where the density partials, and hence det, are ~170x
    different — the multiplier on b1 in dP/dt jumps from ~4e2 to ~7e4), and
    the residual norm jumps 2-3 orders of magnitude and oscillates forever.
    That was diagnosed on the SEGS boiler-internals subsystem: residual fell
    to 1.6e2 by iteration 38, jumped to 1.4e5 by 41, then plateaued at
    ~1.6e5 through iteration 300+ — with EVERY other residual in the network
    still contracting smoothly the whole time, and drum.dP/dt alone carrying
    the entire norm. Setting level_target (0.5 for a half-full drum, the
    usual design point) closes `h` and that same network converges in 89
    iterations to 8.7e-08. Leave at the default None only for
    solve_transient(), where dt drives h forward from h0 instead and no
    steady closure is wanted.

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
        P_target: float | None = None,
        water_out_mdot: float | None = None,
        level_target: float | None = None,
    ):
        if V <= 0.0:
            raise ValueError(f"Drum {name!r}: V must be > 0, got {V}")
        if not (0.0 < level0 < 1.0):
            raise ValueError(f"Drum {name!r}: level0 must be in (0, 1), got {level0}")
        if level_target is not None and not (0.0 < level_target < 1.0):
            raise ValueError(
                f"Drum {name!r}: level_target must be in (0, 1), got {level_target}"
            )
        require_two_phase(fluid, f"Drum {name!r}")

        self.name = name
        self.V = V
        self.has_riser = has_riser
        self.heat_loss = heat_loss
        self.P0 = settings.pressure_to_si(P0)
        self.P_target = None if P_target is None else settings.pressure_to_si(P_target)
        self.water_out_mdot = water_out_mdot
        self.level_target = level_target

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

    @staticmethod
    def _level_from_ph(fluid: "BaseFluid", P: float, h: float, clamp: bool) -> float:
        """Liquid volume fraction at (P, h).

        clamp=True (reporting) pins the mass quality to [0, 1] so a drum
        driven momentarily outside the dome still reports a level in range.
        clamp=False (the level_target residual) must NOT do that: a clamped
        residual is exactly flat in h outside the dome, which puts the
        Jacobian column for this component's `h` back at zero in precisely
        the region the closure exists to keep Newton out of — the unclamped
        form stays monotone in h on both sides of h_f/h_g instead, so a
        Newton step that overshoots the dome still gets pushed back in.
        """
        h_f = fluid.saturated_liquid_enthalpy(P)
        h_g = fluid.saturated_vapor_enthalpy(P)
        x = (h - h_f) / (h_g - h_f) if h_g > h_f else 0.0
        if clamp:
            x = min(max(x, 0.0), 1.0)
        rho = fluid.density_ph(P, h)
        rho_f = fluid.density_ph(P, h_f)
        # level = (1 - x) * v_f / v_avg, with v = 1/rho.
        return (1.0 - x) * rho / rho_f

    @staticmethod
    def _one_sided_eps(
        eps: float, in_dome: Callable[[float], bool], base_in_dome: bool
    ) -> float:
        """`eps`, negated if stepping forward by it would change which side
        of the saturation dome the state sits on.

        A finite-difference step that STRADDLES the dome boundary doesn't
        just lose accuracy — it produces a partial derivative that is a
        blend of two regimes differing by ~10x, and (because the 2x2 solve
        below divides by det = V^2*(rho*drho_dP + drho_dh)) that is enough
        to flip det's SIGN. Measured on the SEGS drum at P=99.98 bar with
        h one J/kg below h_f: the straddling forward step gives
        drho_dh = -2.04e-3 and det = -6.4e-2, against the true one-sided
        liquid-side values drho_dh = -4.3e-4, det = +9.6e-2. dP/dt and dh/dt
        then reverse sign across a 2 J/kg-wide sliver of h, which is a
        discontinuity sitting right where a boiler drum's state naturally
        lives (see _inlet_cp in components/heat_exchanger.py for the same
        class of defect on the heat-exchanger side). Stepping away from the
        boundary instead keeps the derivative one-sided and single-regime;
        backward differencing is just as valid a first-order approximation
        as forward.
        """
        if in_dome(eps) is base_in_dome:
            return eps
        if in_dome(-eps) is base_in_dome:
            return -eps
        return eps  # both directions cross (eps wider than the dome itself)

    def state_derivative(self, state: "NetworkState") -> dict[str, float]:
        P = state.param(f"{self.name}.P")
        h = state.param(f"{self.name}.h")
        fluid = state.fluid_at(self._feed_in_node)
        rho = fluid.density_ph(P, h)

        eps_P = max(abs(P) * 1.0e-6, 1.0)
        eps_h = max(abs(h) * 1.0e-6, 1.0)
        # Keep both perturbations on the same side of the saturation dome as
        # (P, h) itself -- see _one_sided_eps for why a straddling step is
        # sign-flipping rather than merely imprecise.
        try:
            base_in_dome = (
                fluid.saturated_liquid_enthalpy(P) <= h <= fluid.saturated_vapor_enthalpy(P)
            )
            eps_h = self._one_sided_eps(
                eps_h,
                lambda d: (
                    fluid.saturated_liquid_enthalpy(P)
                    <= h + d
                    <= fluid.saturated_vapor_enthalpy(P)
                ),
                base_in_dome,
            )
            eps_P = self._one_sided_eps(
                eps_P,
                lambda d: (
                    fluid.saturated_liquid_enthalpy(P + d)
                    <= h
                    <= fluid.saturated_vapor_enthalpy(P + d)
                ),
                base_in_dome,
            )
        except FluidRangeError:
            # Supercritical / off-chart P has no dome boundary to step away
            # from; the plain forward differences below are already correct.
            pass
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
        h_drum = state.param(f"{self.name}.h")

        P_steam, h_steam = state.node(self._steam_out_node)
        P_water, h_water = state.node(self._water_out_node)
        h_g = fluid.saturated_vapor_enthalpy(P_drum)
        h_f = fluid.saturated_liquid_enthalpy(P_drum)

        out = [
            P_steam - P_drum,
            h_steam - h_g,
            P_water - P_drum,
            h_water - h_f,
        ]
        if self.P_target is not None:
            out.append(P_drum - self.P_target)
        if self.water_out_mdot is not None:
            out.append(state.mdot(self._water_out_node) - self.water_out_mdot)
        if self.level_target is not None:
            out.append(self._level_from_ph(fluid, P_drum, h_drum, clamp=False) - self.level_target)
        return out

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P = state.param(f"{self.name}.P")
        h = state.param(f"{self.name}.h")
        fluid = state.fluid_at(self._feed_in_node)

        level = self._level_from_ph(fluid, P, h, clamp=True)

        return {
            "P [Pa]": P,
            "T_sat [K]": fluid.saturation_temperature(P),
            "level [-]": min(max(level, 0.0), 1.0),
            "mdot_steam [kg/s]": state.mdot(self._steam_out_node),
            "mdot_water [kg/s]": state.mdot(self._water_out_node),
            "V [m^3]": self.V,
        }
