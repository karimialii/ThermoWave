from __future__ import annotations

import math
from typing import TYPE_CHECKING, Callable, Optional

from thermowave.components.base_component import BaseComponent
from thermowave.core.exceptions import FluidRangeError
from thermowave.fluids.two_phase import supports_two_phase

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid

_MIN_C = 1.0e-9  # kg/s * J/(kg K), floor to avoid divide-by-zero if a side has ~0 mdot
_C_SMOOTHING_EPS = 0.002  # kg/s * J/(kg K), smoothing width for _smooth_min's kink
_MIN_MDOT = 1.0e-9  # kg/s, floor for the Q/mdot energy-residual divisions below
_CP_SAT_NUDGE_K = 0.05  # K, how far off the saturation curve _inlet_cp clamps to

_ARRANGEMENTS = ("counterflow", "parallel", "crossflow", "shell_and_tube", "custom")


def _inlet_cp(fluid: "BaseFluid", P: float, T: float) -> float:
    """cp(P, T), clamped at least _CP_SAT_NUDGE_K off the saturation
    temperature for two-phase-capable fluids -- toward whichever side T
    already sits on (liquid if T < T_sat, vapor if T >= T_sat), so a
    legitimately-superheated-vapor inlet elsewhere never gets silently
    treated as subcooled liquid. A real boiler drum's own downcomer
    (Drum.water_out) is ALWAYS exactly saturated liquid by construction --
    not a rare or ill-posed input -- so a two-stream HeatExchanger fed from
    one (as in the SEGS benchmark's evaporator/riser) genuinely, repeatably
    lands its cold inlet right on the saturation boundary at convergence,
    where CoolProp's own phase determination is ambiguous to within its own
    tolerance.

    An earlier version of this only retried ON FAILURE (try cp(T), except
    FluidRangeError: cp(T - nudge), ALWAYS toward liquid) -- that both
    mishandled a genuine vapor-side inlet and, worse, was a DISCONTINUOUS
    function of T (whichever side of CoolProp's own ~1e-4% ambiguity band T
    happens to land on flips which branch runs), which visibly broke
    Newton's own finite-difference Jacobian near this exact boundary: the
    SEGS benchmark's residual norm would drop to within 1% of converging,
    then jump back up by 2 orders of magnitude and oscillate indefinitely --
    the classic signature of a discontinuous residual near the solution,
    not a slow solve. Clamping to a fixed nudge band on the side T is
    ALREADY on (not conditionally, on whether the raw call happened to
    throw) removes that artifact-driven jump; a real jump still exists (cp
    genuinely differs liquid vs. vapor) but it's now fixed to a small,
    whole nudge band on the physically-correct side, not an unpredictable
    ~1e-4%-wide sliver that could land on either.
    """
    if supports_two_phase(fluid):
        try:
            T_sat = fluid.saturation_temperature(P)
        except FluidRangeError:
            # supports_two_phase() is structural (does CoolPropFluid expose these
            # methods at all), not "does this specific backend fluid have a dome" --
            # an incompressible/liquid-only fluid (e.g. INCOMP::TVP1 oil) has no
            # saturation curve to clamp away from at all; nothing to do here.
            T_sat = None
        if T_sat is not None:
            if T_sat - _CP_SAT_NUDGE_K < T < T_sat:
                T = T_sat - _CP_SAT_NUDGE_K
            elif T_sat <= T < T_sat + _CP_SAT_NUDGE_K:
                T = T_sat + _CP_SAT_NUDGE_K
    return fluid.cp(P, T)


class HeatExchanger(BaseComponent):
    """Two-stream heat exchanger: fixed effectiveness, or effectiveness
    derived from geometry (UA, flow arrangement), auto-selected by which one
    you give.

    Four ports: hot_in/hot_out and cold_in/cold_out. Each side's fluid is
    resolved independently via NetworkState.fluid_at() on its own port, so
    two genuinely different fluids (e.g. an oil HTF stream on hot_in/out, the
    network's own working fluid on cold_in/out) work correctly as long as the
    network actually has two fluids to resolve -- see Source's own `fluid`
    parameter for how that second fluid gets introduced in the first place.

    Give at most one of `effectiveness` (fixed-effectiveness mode -- a
    datasheet or an existing exchanger's known performance rating) or `UA`
    (NTU/geometry mode -- effectiveness derived from geometry and a flow
    arrangement, the way a map-based Compressor replaces SimpleCompressor's
    fixed PR/eta). `n_passes`/`arrangement`/`correlation` only apply in UA
    mode; passing any of them alongside `effectiveness` is an error. Giving
    neither also means UA mode, with UA itself left as a free Newton unknown
    (closed externally by a Setpoint/Controller targeting some downstream
    quantity -- typically "T_hot_out [K]" or "T_cold_out [K]" against a
    target approach/pinch temperature, the standard way a real HX's duty
    gets sized to hit a specified terminal temperature difference rather
    than a datasheet UA) -- the same free-parameter pattern SimpleHeatExchanger's
    own Q=None already uses.

    Fixed-effectiveness mode:
        C_hot = mdot_hot * cp_hot, C_cold = mdot_cold * cp_cold (cp evaluated
        at each side's own inlet state), Cmin = min(C_hot, C_cold), and
        Q = effectiveness * Cmin * (T_hot_in - T_cold_in)
    Q is not clamped to be >= 0: a network wired with the "hot" side actually
    colder than the "cold" side just yields Q <= 0 (heat flowing the other
    way), rather than being silently forced to zero.

    UA/NTU mode: single-pass effectiveness comes from the standard
    effectiveness-NTU relations (Cmin = min(mdot*cp) on each side,
    Cr = Cmin/Cmax, NTU = UA/Cmin):
        counterflow: eff = (1-exp(-NTU*(1-Cr))) / (1-Cr*exp(-NTU*(1-Cr)))
                     (eff = NTU/(1+NTU) at Cr == 1)
        parallel:    eff = (1-exp(-NTU*(1+Cr))) / (1+Cr)
        crossflow:   eff = 1 - exp((1/Cr)*NTU**0.22*(exp(-Cr*NTU**0.78)-1))
                     (both fluids unmixed; the standard Incropera correlation)

    For "counterflow"/"parallel"/"crossflow", n_passes > 1 chains that many
    such stages in series along both streams (UA split n_passes ways, one
    set of internal nodes per stream -- same n_elem discretization idea Pipe
    already uses, generalized to two streams), each stage's outlet feeding
    the next stage's inlet. For "counterflow", the two streams are chained
    in opposite physical order (stage 0 is the hot inlet's end and the cold
    outlet's end, the way a real counterflow exchanger's two ends work);
    "parallel"/"crossflow" chain both streams the same direction.

    Be aware of what n_passes does and doesn't buy you for those three
    arrangements: with uniform mdot/cp/UA along the exchanger (the only case
    this component models -- no per-pass property variation), chaining n
    counterflow (or n parallel) stages in series is mathematically
    *identical* to one stage of the same total UA -- subdividing an exact
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
    fluid reversing direction once per shell -- the standard "1-2N" TEMA
    configuration), not a discretization count, and -- unlike the three
    arrangements above -- n_passes genuinely changes the answer here, via a
    real closed-form F-correction-factor effectiveness relation (Bowman,
    Mueller & Nagle 1940, as given in Incropera): one shell pass' own
    effectiveness is
        eff_1 = 2 / (1 + Cr + sqrt(1+Cr^2) * (1+exp(-NTU_1*sqrt(1+Cr^2)))
                                            / (1-exp(-NTU_1*sqrt(1+Cr^2))))
    with NTU_1 = NTU/N (NTU = UA/Cmin for the *whole* exchanger, split evenly
    across the N shells), and N shells combine in series as
        eff = [((1-eff_1*Cr)/(1-eff_1))^N - 1] / [((1-eff_1*Cr)/(1-eff_1))^N - Cr]
    (eff = N*eff_1 / (1+(N-1)*eff_1) at Cr == 1), which reduces to exactly
    eff_1 at N=1 and, as N grows, approaches the true counterflow limit --
    the real physical behavior a genuine multi-pass shell-and-tube exchanger
    has, unlike the other three arrangements' n_passes. No internal
    discretization nodes are created for this arrangement (the tube-side
    reversal is captured analytically by the formula, not modeled as a
    physical node chain) -- see internal_nodes().

    "custom" plugs in a user-supplied effectiveness correlation instead of one
    of the three built-in staged arrangements above -- pass
    `correlation=Callable[[NTU, Cr], float]` (same (NTU, Cr) -> effectiveness
    signature as `_effectiveness()`) for a plate-fin/finned-tube/vendor
    correlation without adding a new named preset here. It chains n_passes
    the same same-direction way "parallel"/"crossflow" do (no arrangement-
    specific directionality is assumed for an arbitrary correlation) and
    splits PR/UA the same per-pass way as the other three staged
    arrangements. `correlation` is required when `arrangement="custom"` and
    must be omitted otherwise.

    Each side's pressure drop is a simple fixed pressure ratio (P_out =
    PR * P_in, same style as SimpleCompressor/SimpleTurbine), not a K-factor
    loss model, in either mode. In UA/staged mode it's split evenly in log
    space across passes (PR_hot**(1/n_passes) per stage) so the overall
    ratio across the whole cascade multiplies back out to PR_hot exactly;
    for "shell_and_tube", applied once across the whole exchanger (there's
    no per-stage node chain to split it across).
    """

    def __init__(
        self,
        name: str,
        PR_hot: float,
        PR_cold: float,
        effectiveness: Optional[float] = None,
        UA: Optional[float] = None,
        n_passes: int = 1,
        arrangement: str = "counterflow",
        correlation: Optional[Callable[[float, float], float]] = None,
        UA_guess: float = 1.0e5,
    ):
        if effectiveness is not None and UA is not None:
            raise ValueError(
                f"HeatExchanger {name!r}: give at most one of effectiveness "
                f"(fixed-effectiveness mode) or UA (NTU/geometry mode) -- got both "
                f"(effectiveness={effectiveness!r}, UA={UA!r})"
            )
        if not (0.0 < PR_hot <= 1.0):
            raise ValueError(f"HeatExchanger {name!r}: PR_hot must be in (0, 1], got {PR_hot}")
        if not (0.0 < PR_cold <= 1.0):
            raise ValueError(f"HeatExchanger {name!r}: PR_cold must be in (0, 1], got {PR_cold}")

        if effectiveness is not None:
            if not (0.0 <= effectiveness <= 1.0):
                raise ValueError(
                    f"HeatExchanger {name!r}: effectiveness must be in [0, 1], "
                    f"got {effectiveness}"
                )
            if n_passes != 1 or arrangement != "counterflow" or correlation is not None:
                raise ValueError(
                    f"HeatExchanger {name!r}: n_passes/arrangement/correlation only apply "
                    f"in UA mode (UA given, or left free); got effectiveness mode instead."
                )
            self._mode = "effectiveness"
        else:
            # UA mode -- UA given (fixed) or left None (a free Newton unknown,
            # closed externally by a Setpoint/Controller targeting some
            # downstream quantity, the same free-parameter pattern
            # SimpleHeatExchanger's own Q=None already uses). Leaving UA free
            # is the genuine equivalent of TESPy's own connection-level
            # pinch/approach-temperature specs (e.g. ttd_l): TESPy's
            # HeatExchanger with no UA/kA/ttd spec at all is *also*
            # under-determined on its own -- it relies on an *additional*
            # target elsewhere (a ttd_l, or another connection's fixed T/P)
            # to close the system, exactly what UA=None + Setpoint(...,
            # target_metric="T_hot_out [K]"/"T_cold_out [K]", ...) gives here.
            if UA is not None and UA <= 0:
                raise ValueError(f"HeatExchanger {name!r}: UA must be > 0, got {UA}")
            if n_passes < 1:
                raise ValueError(f"HeatExchanger {name!r}: n_passes must be >= 1, got {n_passes}")
            if arrangement not in _ARRANGEMENTS:
                raise ValueError(
                    f"HeatExchanger {name!r}: arrangement must be one of "
                    f"{_ARRANGEMENTS}, got {arrangement!r}"
                )
            if arrangement == "custom" and correlation is None:
                raise ValueError(
                    f"HeatExchanger {name!r}: correlation is required when "
                    f"arrangement='custom'"
                )
            if arrangement != "custom" and correlation is not None:
                raise ValueError(
                    f"HeatExchanger {name!r}: correlation is only used when "
                    f"arrangement='custom', got arrangement={arrangement!r}"
                )
            self._mode = "ua"

        self.name = name
        self.effectiveness = effectiveness
        self.UA = UA
        self.UA_guess = UA_guess
        self.PR_hot = PR_hot
        self.PR_cold = PR_cold
        self.n_passes = n_passes
        self.arrangement = arrangement
        self.correlation = correlation
        self._hot_in_node = f"{name}.hot_in"
        self._hot_out_node = f"{name}.hot_out"
        self._cold_in_node = f"{name}.cold_in"
        self._cold_out_node = f"{name}.cold_out"

    def free_parameters(self) -> dict[str, float]:
        if self._mode == "ua" and self.UA is None:
            return {"UA": self.UA_guess}
        return {}

    def _get_UA(self, state: "NetworkState") -> float:
        if self.UA is not None:
            return self.UA
        return state.param(f"{self.name}.UA")

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
        if self._mode == "effectiveness" or self.arrangement == "shell_and_tube":
            # Fixed-effectiveness mode is always a single lumped duty calc
            # (no per-pass chain); shell_and_tube captures its multi-pass
            # behavior analytically -- see this class's own docstring.
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
        # Each pair is warm-started independently (the solver's propagation
        # loop doesn't have both sides' inlet guesses on hand at once to run
        # the real duty calculation), so this is a crude fixed-swing
        # stand-in, not the actual effectiveness-weighted duty: recuperators
        # in this class of machine typically preheat the cold (compressor-
        # discharge) side and cool the hot (turbine-exhaust) side by several
        # hundred K, so a same-order-of-magnitude guess in the right
        # direction is enough to keep downstream map-based components (e.g.
        # a combustor + turbine on the cold side) from starting Newton at a
        # wildly wrong operating point.
        if pair == ("cold_in", "cold_out"):
            return self.PR_cold * P_in, h_in + 4.0e5
        if pair == ("hot_in", "hot_out"):
            return self.PR_hot * P_in, h_in - 3.0e5
        return P_in, h_in

    @staticmethod
    def _smooth_min(a: float, b: float, eps: float = _C_SMOOTHING_EPS) -> float:
        # min(a, b) is exact but non-differentiable at a == b -- a Newton
        # solve landing there (e.g. both sides' warm-started mdot guesses
        # defaulting to the same flat value, which happens easily once mass
        # flow is a free unknown rather than fixed) gets an ill-conditioned
        # or outright singular finite-difference Jacobian right at that
        # point, independent of any variable scaling. This C1-smooth
        # approximation (exact min(a, b) away from the kink, blended over a
        # width of eps right at a == b) is the standard fix for exactly this
        # kind of switching non-smoothness in equation-oriented solvers.
        return 0.5 * (a + b) - 0.5 * ((a - b) ** 2 + eps**2) ** 0.5

    def _effectiveness(self, NTU: float, Cr: float) -> float:
        if self.arrangement == "counterflow":
            if abs(Cr - 1.0) < 1.0e-9:
                return NTU / (1.0 + NTU)
            exp_term = math.exp(-NTU * (1.0 - Cr))
            return (1.0 - exp_term) / (1.0 - Cr * exp_term)
        if self.arrangement == "parallel":
            return (1.0 - math.exp(-NTU * (1.0 + Cr))) / (1.0 + Cr)
        if self.arrangement == "custom":
            assert self.correlation is not None  # enforced in __init__
            return self.correlation(NTU, Cr)
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

    def _duty_fixed_effectiveness(self, state: "NetworkState") -> float:
        P_hot_in, h_hot_in = state.node(self._hot_in_node)
        P_cold_in, h_cold_in = state.node(self._cold_in_node)
        hot_fluid = state.fluid_at(self._hot_in_node)
        cold_fluid = state.fluid_at(self._cold_in_node)
        T_hot_in = hot_fluid.temperature_ph(P_hot_in, h_hot_in)
        T_cold_in = cold_fluid.temperature_ph(P_cold_in, h_cold_in)

        mdot_hot = state.mdot(self._hot_in_node)
        mdot_cold = state.mdot(self._cold_in_node)
        cp_hot = _inlet_cp(hot_fluid, P_hot_in, T_hot_in)
        cp_cold = _inlet_cp(cold_fluid, P_cold_in, T_cold_in)

        C_hot = max(mdot_hot * cp_hot, _MIN_C)
        C_cold = max(mdot_cold * cp_cold, _MIN_C)
        C_min = self._smooth_min(C_hot, C_cold)

        return self.effectiveness * C_min * (T_hot_in - T_cold_in)

    def residuals(self, state: "NetworkState") -> list[float]:
        if self._mode == "effectiveness":
            return self._residuals_fixed_effectiveness(state)
        if self.arrangement == "shell_and_tube":
            return self._residuals_shell_and_tube(state)
        return self._residuals_staged(state)

    def _residuals_fixed_effectiveness(self, state: "NetworkState") -> list[float]:
        Q = self._duty_fixed_effectiveness(state)

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

    def _residuals_shell_and_tube(self, state: "NetworkState") -> list[float]:
        # Single overall duty calc (like fixed-effectiveness mode) rather
        # than the staged per-pass chain below -- the N-shell F-correction
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
        cp_hot = _inlet_cp(hot_fluid, P_hot_in, T_hot_in)
        cp_cold = _inlet_cp(cold_fluid, P_cold_in, T_cold_in)

        C_hot = max(mdot_hot * cp_hot, _MIN_C)
        C_cold = max(mdot_cold * cp_cold, _MIN_C)
        C_min = self._smooth_min(C_hot, C_cold)
        C_max = max(C_hot + C_cold - C_min, C_min)
        Cr = C_min / C_max
        NTU = self._get_UA(state) / C_min
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
        UA_seg = self._get_UA(state) / n
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
            cp_hot = _inlet_cp(hot_fluid, P_hot_a, T_hot_a)
            cp_cold = _inlet_cp(cold_fluid, P_cold_a, T_cold_a)

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
        P_cold_out, h_cold_out = state.node(self._cold_out_node)
        hot_fluid = state.fluid_at(self._hot_in_node)
        cold_fluid = state.fluid_at(self._cold_in_node)
        T_hot_in = hot_fluid.temperature_ph(P_hot_in, h_hot_in)
        T_cold_in = cold_fluid.temperature_ph(P_cold_in, h_cold_in)
        # T_hot_out/T_cold_out use each outlet's own resolved fluid (hot_out's
        # may differ from hot_in's for a composition-changing network, even
        # though this class doesn't itself change composition) rather than
        # assuming hot_fluid/cold_fluid still apply downstream.
        T_hot_out = state.fluid_at(self._hot_out_node).temperature_ph(P_hot_out, h_hot_out)
        T_cold_out = state.fluid_at(self._cold_out_node).temperature_ph(P_cold_out, h_cold_out)

        mdot_hot = state.mdot(self._hot_in_node)
        mdot_cold = state.mdot(self._cold_in_node)
        cp_hot = _inlet_cp(hot_fluid, P_hot_in, T_hot_in)
        cp_cold = _inlet_cp(cold_fluid, P_cold_in, T_cold_in)
        C_min = self._smooth_min(max(mdot_hot * cp_hot, _MIN_C), max(mdot_cold * cp_cold, _MIN_C))
        Q_max = C_min * (T_hot_in - T_cold_in)

        Q = mdot_hot * (h_hot_in - h_hot_out)
        effectiveness = Q / Q_max if abs(Q_max) > _MIN_C else 0.0
        return {
            "power [W]": Q,
            "effectiveness [-]": effectiveness,
            "T_hot_in [K]": T_hot_in,
            "T_cold_in [K]": T_cold_in,
            "T_hot_out [K]": T_hot_out,
            "T_cold_out [K]": T_cold_out,
            "PR_hot [-]": self.PR_hot,
            "PR_cold [-]": self.PR_cold,
        }


class MultiPassHeatExchanger(HeatExchanger):
    """UA/NTU-mode HeatExchanger, kept as its own explicit, always-UA-mode
    class for callers who want that guarantee at the type level. No math is
    duplicated -- see HeatExchanger for the full model.
    """

    def __init__(
        self,
        name: str,
        UA: float,
        PR_hot: float,
        PR_cold: float,
        n_passes: int = 1,
        arrangement: str = "counterflow",
        correlation: Optional[Callable[[float, float], float]] = None,
    ):
        super().__init__(
            name,
            PR_hot,
            PR_cold,
            UA=UA,
            n_passes=n_passes,
            arrangement=arrangement,
            correlation=correlation,
        )
