from __future__ import annotations

import math
from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_MIN_C = 1.0e-9  # kg/s * J/(kg K), floor to avoid divide-by-zero if a side has ~0 mdot
_C_SMOOTHING_EPS = 0.002  # kg/s * J/(kg K), smoothing width for _smooth_min's kink
_MIN_MDOT = 1.0e-9  # kg/s, floor for the Q/mdot energy-residual divisions below

_ARRANGEMENTS = ("counterflow", "parallel", "crossflow", "shell_and_tube")


class MultiPassHeatExchanger(BaseComponent):
    """Two-stream heat exchanger with effectiveness derived from geometry
    (UA, flow arrangement) instead of taken directly as an input — the
    geometry-based counterpart to SimpleHeatExchanger, the way map-based
    Compressor replaced SimpleCompressor's fixed PR/eta.

    Single-pass effectiveness comes from the standard effectiveness-NTU
    relations (Cmin = min(mdot*cp) on each side, Cr = Cmin/Cmax,
    NTU = UA/Cmin):
        counterflow: eff = (1-exp(-NTU*(1-Cr))) / (1-Cr*exp(-NTU*(1-Cr)))
                     (eff = NTU/(1+NTU) at Cr == 1)
        parallel:    eff = (1-exp(-NTU*(1+Cr))) / (1+Cr)
        crossflow:   eff = 1 - exp((1/Cr)*NTU**0.22*(exp(-Cr*NTU**0.78)-1))
                     (both fluids unmixed; the standard Incropera correlation)

    For "counterflow"/"parallel"/"crossflow", n_passes > 1 chains that many
    such stages in series along both streams (UA split n_passes ways, one
    set of internal nodes per stream — same n_elem discretization idea Pipe
    already uses, generalized to two streams), each stage's outlet feeding
    the next stage's inlet. For "counterflow", the two streams are chained
    in opposite physical order (stage 0 is the hot inlet's end and the cold
    outlet's end, the way a real counterflow exchanger's two ends work);
    "parallel"/"crossflow" chain both streams the same direction.

    Be aware of what n_passes does and doesn't buy you for those three
    arrangements: with uniform mdot/cp/UA along the exchanger (the only case
    this component models — no per-pass property variation), chaining n
    counterflow (or n parallel) stages in series is mathematically
    *identical* to one stage of the same total UA — subdividing an exact
    closed-form solution along a fixed flow direction doesn't change it, so
    n_passes=1..N give the same effectiveness for those two arrangements.
    It's not a no-op for "crossflow": the Incropera correlation above is
    itself only an approximation of true (2D) crossflow heat transfer, and
    chaining same-direction 1D slices of it is a further approximation on
    top of that one, so its result under subdivision isn't guaranteed to
    move monotonically toward counterflow's performance. Treat n_passes
    mainly as node/reporting granularity for these three arrangements
    (useful groundwork for a future per-pass-property extension) rather than
    a guaranteed performance lever; n_passes=1 is the well-defined, exact
    case for every arrangement.

    "shell_and_tube" is the genuine reversing-header multi-pass case (real
    shell-and-tube exchanger hardware): n_passes here means the number of
    *shell* passes N (each shell pass internally has 2 tube passes, tube
    fluid reversing direction once per shell — the standard "1-2N" TEMA
    configuration), not a discretization count, and — unlike the three
    arrangements above — n_passes genuinely changes the answer here, via a
    real closed-form F-correction-factor effectiveness relation (Bowman,
    Mueller & Nagle 1940, as given in Incropera): one shell pass' own
    effectiveness is
        eff_1 = 2 / (1 + Cr + sqrt(1+Cr^2) * (1+exp(-NTU_1*sqrt(1+Cr^2)))
                                            / (1-exp(-NTU_1*sqrt(1+Cr^2))))
    with NTU_1 = NTU/N (NTU = UA/Cmin for the *whole* exchanger, split evenly
    across the N shells), and N shells combine in series as
        eff = [((1-eff_1*Cr)/(1-eff_1))^N - 1] / [((1-eff_1*Cr)/(1-eff_1))^N - Cr]
    (eff = N*eff_1 / (1+(N-1)*eff_1) at Cr == 1), which reduces to exactly
    eff_1 at N=1 and, as N grows, approaches the true counterflow limit —
    the real physical behavior a genuine multi-pass shell-and-tube exchanger
    has, unlike the other three arrangements' n_passes. No internal
    discretization nodes are created for this arrangement (the tube-side
    reversal is captured analytically by the formula, not modeled as a
    physical node chain) — see internal_nodes().

    Each side's pressure drop is still a fixed ratio (PR_hot/PR_cold, same
    style as SimpleHeatExchanger) — for "counterflow"/"parallel"/
    "crossflow", split evenly in log space across passes
    (PR_hot**(1/n_passes) per stage) so the overall ratio across the whole
    cascade multiplies back out to PR_hot exactly; for "shell_and_tube",
    applied once across the whole exchanger (there's no per-stage node chain
    to split it across).
    """

    def __init__(
        self,
        name: str,
        UA: float,
        PR_hot: float,
        PR_cold: float,
        n_passes: int = 1,
        arrangement: str = "counterflow",
    ):
        if UA <= 0:
            raise ValueError(f"MultiPassHeatExchanger {name!r}: UA must be > 0, got {UA}")
        if PR_hot <= 0:
            raise ValueError(f"MultiPassHeatExchanger {name!r}: PR_hot must be > 0, got {PR_hot}")
        if PR_cold <= 0:
            raise ValueError(
                f"MultiPassHeatExchanger {name!r}: PR_cold must be > 0, got {PR_cold}"
            )
        if n_passes < 1:
            raise ValueError(
                f"MultiPassHeatExchanger {name!r}: n_passes must be >= 1, got {n_passes}"
            )
        if arrangement not in _ARRANGEMENTS:
            raise ValueError(
                f"MultiPassHeatExchanger {name!r}: arrangement must be one of "
                f"{_ARRANGEMENTS}, got {arrangement!r}"
            )
        self.name = name
        self.UA = UA
        self.PR_hot = PR_hot
        self.PR_cold = PR_cold
        self.n_passes = n_passes
        self.arrangement = arrangement
        self._hot_in_node = f"{name}.hot_in"
        self._hot_out_node = f"{name}.hot_out"
        self._cold_in_node = f"{name}.cold_in"
        self._cold_out_node = f"{name}.cold_out"

    def ports(self) -> dict[str, str]:
        return {
            "hot_in": self._hot_in_node,
            "hot_out": self._hot_out_node,
            "cold_in": self._cold_in_node,
            "cold_out": self._cold_out_node,
        }

    def report_category(self) -> str:
        return "heat_exchanger"

    def internal_nodes(self) -> list[str]:
        if self.arrangement == "shell_and_tube":
            # The tube-side reversal per shell pass is captured analytically
            # by _effectiveness_shell_and_tube() below, not modeled as a
            # physical node chain -- see this class's own docstring.
            return []
        return [f"{self.name}__hot_mid{i}" for i in range(1, self.n_passes)] + [
            f"{self.name}__cold_mid{i}" for i in range(1, self.n_passes)
        ]

    def _hot_chain(self) -> list[str]:
        mids = [f"{self.name}__hot_mid{i}" for i in range(1, self.n_passes)]
        return [self._hot_in_node, *mids, self._hot_out_node]

    def _cold_chain(self) -> list[str]:
        mids = [f"{self.name}__cold_mid{i}" for i in range(1, self.n_passes)]
        return [self._cold_in_node, *mids, self._cold_out_node]

    def warm_start_pairs(self) -> list[tuple[str, str]]:
        return [("hot_in", "hot_out"), ("cold_in", "cold_out")]

    def guess_outlet_for_pair(
        self, pair: tuple[str, str], P_in: float, h_in: float, mdot: float
    ) -> tuple[float, float]:
        # Same crude fixed-swing rationale as SimpleHeatExchanger's own
        # override — see its docstring: each pair is warm-started
        # independently, without both sides' inlet guesses on hand at once
        # to run the real per-pass duty calculation.
        if pair == ("cold_in", "cold_out"):
            return self.PR_cold * P_in, h_in + 4.0e5
        if pair == ("hot_in", "hot_out"):
            return self.PR_hot * P_in, h_in - 3.0e5
        return P_in, h_in

    @staticmethod
    def _smooth_min(a: float, b: float, eps: float = _C_SMOOTHING_EPS) -> float:
        # See SimpleHeatExchanger._smooth_min's docstring for why this is
        # needed instead of the exact min(a, b).
        return 0.5 * (a + b) - 0.5 * ((a - b) ** 2 + eps**2) ** 0.5

    def _effectiveness(self, NTU: float, Cr: float) -> float:
        if self.arrangement == "counterflow":
            if abs(Cr - 1.0) < 1.0e-9:
                return NTU / (1.0 + NTU)
            exp_term = math.exp(-NTU * (1.0 - Cr))
            return (1.0 - exp_term) / (1.0 - Cr * exp_term)
        if self.arrangement == "parallel":
            return (1.0 - math.exp(-NTU * (1.0 + Cr))) / (1.0 + Cr)
        # crossflow, both fluids unmixed (Incropera's correlation)
        if Cr < 1.0e-9:
            return 1.0 - math.exp(-NTU)
        return 1.0 - math.exp((1.0 / Cr) * NTU**0.22 * (math.exp(-Cr * NTU**0.78) - 1.0))

    @staticmethod
    def _effectiveness_one_shell_pass(NTU1: float, Cr: float) -> float:
        """eff_1: one shell pass' own effectiveness (Bowman/Mueller/Nagle,
        as given in Incropera) -- see this class's own docstring."""
        root = math.sqrt(1.0 + Cr**2)
        exp_term = math.exp(-NTU1 * root)
        return 2.0 / (1.0 + Cr + root * (1.0 + exp_term) / (1.0 - exp_term))

    @classmethod
    def _effectiveness_shell_and_tube(cls, NTU: float, Cr: float, n_shells: int) -> float:
        """N shell passes in series, each with NTU_1 = NTU/N -- see this
        class's own docstring for the closed-form combination."""
        eff1 = cls._effectiveness_one_shell_pass(NTU / n_shells, Cr)
        if abs(Cr - 1.0) < 1.0e-9:
            return n_shells * eff1 / (1.0 + (n_shells - 1.0) * eff1)
        ratio = ((1.0 - eff1 * Cr) / (1.0 - eff1)) ** n_shells
        return (ratio - 1.0) / (ratio - Cr)

    def _pass_nodes(self, i: int, hot_chain: list[str], cold_chain: list[str]) -> tuple:
        n = self.n_passes
        hot_a, hot_b = hot_chain[i], hot_chain[i + 1]
        if self.arrangement == "counterflow":
            cold_a, cold_b = cold_chain[n - i - 1], cold_chain[n - i]
        else:
            cold_a, cold_b = cold_chain[i], cold_chain[i + 1]
        return hot_a, hot_b, cold_a, cold_b

    def residuals(self, state: "NetworkState") -> list[float]:
        if self.arrangement == "shell_and_tube":
            return self._residuals_shell_and_tube(state)
        return self._residuals_staged(state)

    def _residuals_shell_and_tube(self, state: "NetworkState") -> list[float]:
        # Single overall duty calc (like SimpleHeatExchanger) rather than
        # the staged per-pass chain below -- the N-shell F-correction
        # formula already accounts for the whole multi-pass exchanger
        # analytically, so there's no separate per-stage node chain to
        # solve through (see internal_nodes() and this class's own
        # docstring).
        P_hot_in, h_hot_in = state.node(self._hot_in_node)
        P_cold_in, h_cold_in = state.node(self._cold_in_node)
        hot_fluid = state.fluid_at(self._hot_in_node)
        cold_fluid = state.fluid_at(self._cold_in_node)
        T_hot_in = hot_fluid.temperature_ph(P_hot_in, h_hot_in)
        T_cold_in = cold_fluid.temperature_ph(P_cold_in, h_cold_in)

        mdot_hot = state.mdot(self._hot_in_node)
        mdot_cold = state.mdot(self._cold_in_node)
        cp_hot = hot_fluid.cp(P_hot_in, T_hot_in)
        cp_cold = cold_fluid.cp(P_cold_in, T_cold_in)

        C_hot = max(mdot_hot * cp_hot, _MIN_C)
        C_cold = max(mdot_cold * cp_cold, _MIN_C)
        C_min = self._smooth_min(C_hot, C_cold)
        C_max = max(C_hot + C_cold - C_min, C_min)
        Cr = C_min / C_max
        NTU = self.UA / C_min
        eff = self._effectiveness_shell_and_tube(NTU, Cr, self.n_passes)
        Q = eff * C_min * (T_hot_in - T_cold_in)

        P_hot_out, h_hot_out = state.node(self._hot_out_node)
        P_cold_out, h_cold_out = state.node(self._cold_out_node)

        hot_momentum_residual = P_hot_out - self.PR_hot * P_hot_in
        hot_energy_residual = h_hot_out - (h_hot_in - Q / max(mdot_hot, _MIN_MDOT))
        hot_mass_residual = state.mdot(self._hot_out_node) - mdot_hot

        cold_momentum_residual = P_cold_out - self.PR_cold * P_cold_in
        cold_energy_residual = h_cold_out - (h_cold_in + Q / max(mdot_cold, _MIN_MDOT))
        cold_mass_residual = state.mdot(self._cold_out_node) - mdot_cold

        return [
            hot_momentum_residual,
            hot_energy_residual,
            hot_mass_residual,
            cold_momentum_residual,
            cold_energy_residual,
            cold_mass_residual,
        ]

    def _residuals_staged(self, state: "NetworkState") -> list[float]:
        hot_chain = self._hot_chain()
        cold_chain = self._cold_chain()
        n = self.n_passes
        mdot_hot = state.mdot(self._hot_in_node)
        mdot_cold = state.mdot(self._cold_in_node)
        UA_seg = self.UA / n
        PR_hot_seg = self.PR_hot ** (1.0 / n)
        PR_cold_seg = self.PR_cold ** (1.0 / n)
        # One fluid per stream for the whole exchanger (its own inlet
        # port's, not each per-pass internal node's) -- internal nodes never
        # get their own node_fluid entry (only port nodes do, see
        # Network._resolve_node_fluid()), and neither stream changes
        # composition pass-to-pass anyway.
        hot_fluid = state.fluid_at(self._hot_in_node)
        cold_fluid = state.fluid_at(self._cold_in_node)

        out: list[float] = []
        for i in range(n):
            hot_a, hot_b, cold_a, cold_b = self._pass_nodes(i, hot_chain, cold_chain)
            P_hot_a, h_hot_a = state.node(hot_a)
            P_cold_a, h_cold_a = state.node(cold_a)
            T_hot_a = hot_fluid.temperature_ph(P_hot_a, h_hot_a)
            T_cold_a = cold_fluid.temperature_ph(P_cold_a, h_cold_a)
            cp_hot = hot_fluid.cp(P_hot_a, T_hot_a)
            cp_cold = cold_fluid.cp(P_cold_a, T_cold_a)

            C_hot = max(mdot_hot * cp_hot, _MIN_C)
            C_cold = max(mdot_cold * cp_cold, _MIN_C)
            C_min = self._smooth_min(C_hot, C_cold)
            C_max = max(C_hot + C_cold - C_min, C_min)
            Cr = C_min / C_max
            NTU = UA_seg / C_min
            eff = self._effectiveness(NTU, Cr)
            Q = eff * C_min * (T_hot_a - T_cold_a)

            P_hot_b, h_hot_b = state.node(hot_b)
            P_cold_b, h_cold_b = state.node(cold_b)
            out.append(P_hot_b - PR_hot_seg * P_hot_a)
            out.append(h_hot_b - (h_hot_a - Q / max(mdot_hot, _MIN_MDOT)))
            out.append(P_cold_b - PR_cold_seg * P_cold_a)
            out.append(h_cold_b - (h_cold_a + Q / max(mdot_cold, _MIN_MDOT)))

        out.append(state.mdot(self._hot_out_node) - mdot_hot)
        out.append(state.mdot(self._cold_out_node) - mdot_cold)
        return out

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_hot_in, h_hot_in = state.node(self._hot_in_node)
        P_hot_out, h_hot_out = state.node(self._hot_out_node)
        P_cold_in, h_cold_in = state.node(self._cold_in_node)
        hot_fluid = state.fluid_at(self._hot_in_node)
        cold_fluid = state.fluid_at(self._cold_in_node)
        T_hot_in = hot_fluid.temperature_ph(P_hot_in, h_hot_in)
        T_cold_in = cold_fluid.temperature_ph(P_cold_in, h_cold_in)

        mdot_hot = state.mdot(self._hot_in_node)
        mdot_cold = state.mdot(self._cold_in_node)
        cp_hot = hot_fluid.cp(P_hot_in, T_hot_in)
        cp_cold = cold_fluid.cp(P_cold_in, T_cold_in)
        C_min = self._smooth_min(max(mdot_hot * cp_hot, _MIN_C), max(mdot_cold * cp_cold, _MIN_C))
        Q_max = C_min * (T_hot_in - T_cold_in)

        Q = mdot_hot * (h_hot_in - h_hot_out)
        return {
            "power [W]": Q,
            "effectiveness [-]": Q / Q_max if abs(Q_max) > _MIN_C else 0.0,
            "T_hot_in [K]": T_hot_in,
            "T_cold_in [K]": T_cold_in,
            "PR_hot [-]": self.PR_hot,
            "PR_cold [-]": self.PR_cold,
        }
