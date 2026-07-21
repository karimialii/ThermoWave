from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_MIN_C = 1.0e-9  # kg/s * J/(kg K), floor to avoid divide-by-zero if a side has ~0 mdot
_C_SMOOTHING_EPS = 0.002  # kg/s * J/(kg K), smoothing width for _smooth_min's kink
_MIN_MDOT = 1.0e-9  # kg/s, floor for the Q/mdot energy-residual divisions below —
# unlike C_hot/C_cold (which only feed the effectiveness-NTU duty calc, and are
# already floored by _MIN_C), these divide Q directly, so a Newton iterate that
# clamps a free mdot near zero (Solver.MDOT_MIN) would otherwise blow up here.


class SimpleHeatExchanger(BaseComponent):
    """Two-stream heat exchanger, 0D model with a fixed effectiveness.

    Four ports: hot_in/hot_out and cold_in/cold_out, both streams sharing the
    Network's single fluid model (no distinct hot/cold fluids yet — same
    limitation as every other component here).

    Heat duty uses the effectiveness directly as an input (rather than being
    derived from UA/NTU and a flow arrangement — a fixed rating, the way a
    datasheet or an existing exchanger's known performance would be given):
        C_hot = mdot_hot * cp_hot, C_cold = mdot_cold * cp_cold (cp evaluated
        at each side's own inlet state), Cmin = min(C_hot, C_cold), and
        Q = effectiveness * Cmin * (T_hot_in - T_cold_in)
    Q is not clamped to be >= 0: a network wired with the "hot" side actually
    colder than the "cold" side just yields Q <= 0 (heat flowing the other
    way), rather than being silently forced to zero.

    Each stream's pressure drop is a simple fixed pressure ratio (P_out =
    PR * P_in, same style as SimpleCompressor/SimpleTurbine), not a K-factor
    loss model — geometry-based resistance (and, if U/area/other physical
    parameters are known, computing UA/NTU from them instead of taking
    effectiveness directly) can be layered in later without changing this
    component's interface.
    """

    def __init__(
        self,
        name: str,
        effectiveness: float,
        PR_hot: float,
        PR_cold: float,
    ):
        if not (0.0 <= effectiveness <= 1.0):
            raise ValueError(
                f"SimpleHeatExchanger {name!r}: effectiveness must be in [0, 1], "
                f"got {effectiveness}"
            )
        if PR_hot <= 0:
            raise ValueError(f"SimpleHeatExchanger {name!r}: PR_hot must be > 0, got {PR_hot}")
        if PR_cold <= 0:
            raise ValueError(f"SimpleHeatExchanger {name!r}: PR_cold must be > 0, got {PR_cold}")
        self.name = name
        self.effectiveness = effectiveness
        self.PR_hot = PR_hot
        self.PR_cold = PR_cold
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

    def warm_start_pairs(self) -> list[tuple[str, str]]:
        return [("hot_in", "hot_out"), ("cold_in", "cold_out")]

    def guess_outlet_for_pair(
        self, pair: tuple[str, str], P_in: float, h_in: float, mdot: float
    ) -> tuple[float, float]:
        # Each pair is warm-started independently (the solver's propagation
        # loop doesn't have both sides' inlet guesses on hand at once to run
        # the real _duty() calculation), so this is a crude fixed-swing
        # stand-in, not the actual effectiveness-weighted duty: recuperators
        # in this class of machine (see 26_t100_power_profile.py) typically
        # preheat the cold (compressor-discharge) side and cool the hot
        # (turbine-exhaust) side by several hundred K, so a same-order-of-
        # magnitude guess in the right direction is enough to keep
        # downstream map-based components (e.g. a combustor + turbine on the
        # cold side) from starting Newton at a wildly wrong operating point.
        if pair == ("cold_in", "cold_out"):
            return self.PR_cold * P_in, h_in + 4.0e5
        if pair == ("hot_in", "hot_out"):
            return self.PR_hot * P_in, h_in - 3.0e5
        return P_in, h_in

    @staticmethod
    def _smooth_min(a: float, b: float, eps: float = _C_SMOOTHING_EPS) -> float:
        # min(a, b) is exact but non-differentiable at a == b — a Newton
        # solve landing there (e.g. both sides' warm-started mdot guesses
        # defaulting to the same flat value, which happens easily once mass
        # flow is a free unknown rather than fixed) gets an ill-conditioned
        # or outright singular finite-difference Jacobian right at that
        # point, independent of any variable scaling. This C1-smooth
        # approximation (exact min(a, b) away from the kink, blended over a
        # width of eps right at a == b) is the standard fix for exactly this
        # kind of switching non-smoothness in equation-oriented solvers.
        return 0.5 * (a + b) - 0.5 * ((a - b) ** 2 + eps**2) ** 0.5

    def _duty(self, state: "NetworkState") -> float:
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

        return self.effectiveness * C_min * (T_hot_in - T_cold_in)

    def residuals(self, state: "NetworkState") -> list[float]:
        Q = self._duty(state)

        P_hot_in, h_hot_in = state.node(self._hot_in_node)
        P_hot_out, h_hot_out = state.node(self._hot_out_node)
        P_cold_in, h_cold_in = state.node(self._cold_in_node)
        P_cold_out, h_cold_out = state.node(self._cold_out_node)

        mdot_hot = state.mdot(self._hot_in_node)
        mdot_cold = state.mdot(self._cold_in_node)

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

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_hot_in, h_hot_in = state.node(self._hot_in_node)
        P_cold_in, h_cold_in = state.node(self._cold_in_node)
        return {
            "power [W]": self._duty(state),
            "T_hot_in [K]": state.fluid_at(self._hot_in_node).temperature_ph(P_hot_in, h_hot_in),
            "T_cold_in [K]": state.fluid_at(self._cold_in_node).temperature_ph(P_cold_in, h_cold_in),
            "PR_hot [-]": self.PR_hot,
            "PR_cold [-]": self.PR_cold,
        }
