from __future__ import annotations

import math
from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class Nozzle(BaseComponent):
    """Converging (or, with D_exit, converging-diverging) nozzle: converts a
    pressure/enthalpy drop into velocity — a Turbine without a shaft (no free
    speed parameter, no mechanical power; the enthalpy drop becomes kinetic
    energy instead of shaft work).

    Isentropic expansion from the inlet's (P_in, h_in) — treated as
    stagnation conditions, negligible inlet velocity, the same convention
    map-based Compressor/Turbine already use for their own P03/T03 — down to
    the outlet node's static pressure P_out (fixed by whatever's downstream,
    e.g. a Sink): T_out_isentropic = T_in * (P_out/P_in)**((gamma-1)/gamma),
    dh_ideal = h_in - h(P_out, T_out_isentropic), dh_actual = eta * dh_ideal
    (eta: velocity coefficient / nozzle efficiency, same role as Turbine's
    eta_s), V_exit = sqrt(2 * dh_actual). gamma: give it directly, or leave
    it None (the default) to derive it from the network's own fluid model
    (BaseFluid.gamma()), same convention as Compressor/Turbine/
    SimpleCompressor/SimpleTurbine.

    Choked flow: if P_out/P_in falls below the critical ratio (2/(gamma+1))
    ** (gamma/(gamma-1)), the throat is sonic and dh_ideal/mdot/V_exit are
    capped at their critical-pressure-ratio values instead of continuing to
    grow as P_out drops further — a real, useful limit (predicting the
    maximum mass flow a converging nozzle of a given throat area can pass),
    not just a numerical safeguard.

    Diverging section (D_exit): leave it None (default) for a plain
    converging nozzle, where once choked the true exit *static* pressure is
    just the critical pressure at the throat — the actual external
    supersonic expansion of an under/over-expanded jet down to ambient isn't
    modeled (that's genuinely complex shock/expansion gas dynamics, out of
    scope for a lumped 1D component), so the outlet node's reported (P, h)
    stops representing a single fully self-consistent static state below the
    critical pressure ratio, the same kind of explicit, documented
    simplification CheckValve makes for its own reverse-flow regime.

    Give D_exit (> D, the throat diameter) for a converging-diverging (de
    Laval) nozzle instead: once choked at the throat, the flow continues to
    expand isentropically through the diverging section, reaching whatever
    supersonic Mach number the geometric area ratio A_exit/A_throat implies
    (the standard 1D isentropic area-Mach relation, solved for its
    supersonic root — see _mach_from_area_ratio()). This is the nozzle's
    *design* exit condition: the area ratio alone determines the fully-
    expanded exit static pressure (reported as "P_exit_ideal [Pa]"),
    assuming shock-free isentropic flow all the way to the exit plane. Real
    off-design behavior — a normal shock inside the diverging section for a
    higher back pressure than that design value, or external oblique
    shocks/expansion fans for a physical exit pressure mismatch at very low
    back pressure — isn't modeled (same class of scope boundary as the
    plain converging case above); compare the reported "P_exit_ideal [Pa]"
    against the network's own downstream pressure to judge how far off
    design the current operating point is. mdot is unaffected either way:
    it's set entirely by continuity at the (always sonic-when-choked) throat
    once choked, not by whatever the diverging section does afterward — a
    real nozzle's mass flow can't depend on what's downstream of where it's
    already sonic. When *not* choked, D_exit is inert (the whole nozzle just
    behaves as the plain converging case — a subsonic converging-diverging
    nozzle's diverging section would act as a diffuser, a different flow
    regime this class doesn't attempt to model).

    mdot is not an independent input here (unlike Pipe/Valve, where mdot is
    given and dp is solved for): a nozzle's throat area *determines* mass
    flow from the upstream/downstream pressures, so residuals() instead
    solves for mdot from continuity at the throat (density there × area ×
    velocity), the reverse direction from Pipe's own dp = f(mdot).
    """

    def __init__(
        self,
        name: str,
        D: float,
        eta: float = 0.98,
        gamma: float | None = None,
        D_exit: float | None = None,
    ):
        if not (0.0 < eta <= 1.0):
            raise ValueError(f"Nozzle {name!r}: eta must be in (0, 1], got {eta}")
        if gamma is not None and gamma <= 1.0:
            raise ValueError(f"Nozzle {name!r}: gamma must be > 1, got {gamma}")
        if D_exit is not None and D_exit <= D:
            raise ValueError(
                f"Nozzle {name!r}: D_exit must be > D (the throat diameter) to form a "
                f"diverging section, got D_exit={D_exit} <= D={D}"
            )
        self.name = name
        self.D = D
        self.eta = eta
        self.gamma = gamma
        self.D_exit = D_exit
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"
        self._area = math.pi * D**2 / 4
        self._area_exit = math.pi * D_exit**2 / 4 if D_exit is not None else None

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def report_category(self) -> str:
        return "nozzle"

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        return 0.5 * P_in, h_in

    def _gamma(self, state: "NetworkState", P_in: float, T_in: float) -> float:
        if self.gamma is not None:
            return self.gamma
        return state.fluid_at(self._inlet_node).gamma(P_in, T_in)

    @staticmethod
    def _mach_from_area_ratio(gamma: float, area_ratio: float) -> float:
        """Supersonic root M > 1 of the isentropic area-Mach relation
        A/A* = (1/M) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * M**2))
               ** ((gamma+1) / (2*(gamma-1)))
        for a converging-diverging nozzle's given exit/throat area ratio.
        Solved by bisection: the supersonic branch (M > 1) is monotonically
        increasing in M for a fixed area_ratio (unlike the subsonic branch
        M < 1, which isn't — but a choked throat's downstream flow is always
        on the supersonic branch), and this only needs solving once per
        residual evaluation (not on any hot per-Jacobian-column path), so a
        closed-form/Newton shortcut isn't worth the added complexity.
        """

        def area_ratio_of(M: float) -> float:
            return (1.0 / M) * (
                (2.0 / (gamma + 1.0)) * (1.0 + (gamma - 1.0) / 2.0 * M**2)
            ) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))

        lo, hi = 1.0 + 1.0e-9, 2.0
        while area_ratio_of(hi) < area_ratio and hi < 1.0e6:
            hi *= 2.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if area_ratio_of(mid) < area_ratio:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def _flow(self, state: "NetworkState") -> dict[str, float]:
        """Throat state (always computed — this is what sets mdot via
        continuity) plus, when D_exit is given and the throat is choked,
        the diverging-section exit-plane state too (see this class's own
        docstring). Returns a dict rather than a tuple since the shape
        varies (exit-plane keys only appear when applicable) and both
        residuals() and report_metrics() need different subsets of it.
        """
        P_in, h_in = state.node(self._inlet_node)
        P_out, _ = state.node(self._outlet_node)
        fluid = state.fluid_at(self._inlet_node)

        T_in = fluid.temperature_ph(P_in, h_in)
        gamma = self._gamma(state, P_in, T_in)
        PR = P_out / P_in
        PR_crit = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
        choked = PR < PR_crit
        PR_throat = max(PR, PR_crit)

        T_throat_isentropic = T_in * PR_throat ** ((gamma - 1.0) / gamma)
        h_throat_isentropic = fluid.enthalpy_pt(PR_throat * P_in, T_throat_isentropic)
        dh_throat_ideal = h_in - h_throat_isentropic
        dh_throat_actual = self.eta * dh_throat_ideal
        V_throat = math.sqrt(max(2.0 * dh_throat_actual, 0.0))
        h_throat_actual = h_in - dh_throat_actual
        P_throat = PR_throat * P_in
        rho_throat = fluid.density_ph(P_throat, h_throat_actual)
        mdot_computed = rho_throat * self._area * V_throat
        a_throat = math.sqrt(gamma * P_throat / rho_throat)

        result = {
            "mdot": mdot_computed,
            "h_out": h_throat_actual,
            "PR": PR,
            "choked": choked,
            "V_throat": V_throat,
            "Mach_throat": V_throat / a_throat if a_throat > 0.0 else 0.0,
        }

        if self._area_exit is not None and choked:
            area_ratio = self._area_exit / self._area
            M_exit = self._mach_from_area_ratio(gamma, area_ratio)
            T_exit_isentropic = T_in / (1.0 + (gamma - 1.0) / 2.0 * M_exit**2)
            P_exit_isentropic = P_in * (T_exit_isentropic / T_in) ** (gamma / (gamma - 1.0))
            h_exit_isentropic = fluid.enthalpy_pt(P_exit_isentropic, T_exit_isentropic)
            dh_exit_ideal = h_in - h_exit_isentropic
            dh_exit_actual = self.eta * dh_exit_ideal
            V_exit = math.sqrt(max(2.0 * dh_exit_actual, 0.0))
            h_exit_actual = h_in - dh_exit_actual
            rho_exit = fluid.density_ph(P_exit_isentropic, h_exit_actual)
            a_exit = math.sqrt(gamma * P_exit_isentropic / rho_exit)

            result["h_out"] = h_exit_actual
            result["V_exit"] = V_exit
            result["Mach_exit"] = V_exit / a_exit if a_exit > 0.0 else 0.0
            result["P_exit_ideal"] = P_exit_isentropic

        return result

    def residuals(self, state: "NetworkState") -> list[float]:
        flow = self._flow(state)
        _, h_out = state.node(self._outlet_node)
        mdot_in = state.mdot(self._inlet_node)

        mass_residual = mdot_in - flow["mdot"]
        energy_residual = h_out - flow["h_out"]
        mdot_out_residual = state.mdot(self._outlet_node) - mdot_in
        return [mass_residual, energy_residual, mdot_out_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        flow = self._flow(state)
        metrics = {
            "V [m/s]": flow["V_throat"],
            "Mach [-]": flow["Mach_throat"],
            "mdot [kg/s]": flow["mdot"],
            "PR [-]": flow["PR"],
            "choked [-]": 1.0 if flow["choked"] else 0.0,
        }
        if "V_exit" in flow:
            metrics["V_exit [m/s]"] = flow["V_exit"]
            metrics["Mach_exit [-]"] = flow["Mach_exit"]
            metrics["P_exit_ideal [Pa]"] = flow["P_exit_ideal"]
        return metrics
