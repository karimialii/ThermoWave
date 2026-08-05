from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid


class Recycle(BaseComponent):
    """Mass-flow anchor for a genuinely closed loop -- no Source/Sink needed.

    Every flow component's own mass residual (mdot_out - mdot_in = 0) already
    chains mdot all the way around a loop; what a loop with no Source in it
    is missing isn't a *relationship* between its mdot unknowns (those are
    already tied together), it's an anchor for their shared *magnitude* --
    conservation alone only pins mdot up to an arbitrary constant. Likewise,
    P and h at every node in a real closed cycle (ORC, Rankine, Brayton) are
    already fully determined by the loop's own component physics -- a
    SimpleCondenser computes its own outlet (P, h) from subcool/quality, a
    Pump raises P by a given ratio, etc. -- so nothing needs a Source's other
    job (fixing (P, h) as a hard boundary condition) either.

    Recycle supplies exactly the piece that's missing and nothing else: drop
    it anywhere in a loop (typically condenser outlet -> pump inlet) to pin
    the loop's overall mass-flow scale, with every other node's state left
    entirely to the rest of the network's own physics:

        net.connect(condenser, "out", recycle, "in")
        net.connect(recycle, "out", pump, "in")

    Ports: in, out.

    Fixes: mdot at out only -- never (P, h). This is deliberately the
    opposite half of Source's job: Source fixes (P, h) always and mdot
    optionally, Recycle fixes mdot always and never touches (P, h).

    Residuals: exactly 2 -- P_out - P_in = 0, h_out - h_in = 0, an ideal
    lossless pass-through. NOT the 3-residual pattern most flow components
    use: there's deliberately no mass-conservation residual here, since
    fixed_node_mdot() already fixes this component's own outlet mdot
    directly -- adding a mass residual on top of that would try to pin the
    same already-fixed unknown with a second equation.

    Composition (EGR-style loops): pass fluid_guess=<a BaseFluid> to also
    tear a fluid-*composition* cycle open here, not just the flow cycle.
    Network._resolve_node_fluid() only ever forward-propagates composition
    from the network's own fixed (P, h) boundaries (Source, etc.) -- a loop
    closed by Recycle fixes none of those, so without a guess a node whose
    composition depends on itself through the loop (e.g. a Combustor whose
    inlet is fed, via a Junction, by its own recirculated exhaust) would
    silently fall back to the network's plain default fluid instead of ever
    resolving. fluid_guess seeds this component's outlet with a starting
    composition; Network.solve() then re-solves and refines that guess
    against what actually arrives back at this component's inlet (direct-
    substitution / "tear stream" iteration, the standard technique process
    simulators use for recycle convergence) until it self-consistently
    converges (see Network.solve()'s recycle_fluid_tol/recycle_fluid_max_iter).
    Leave fluid_guess unset (the default) for a loop with no composition-
    changing component in it (e.g. a plain Rankine/Brayton cycle) -- that
    case needs none of this and behaves exactly as before.
    """

    def __init__(
        self,
        name: str,
        mdot: float,
        fluid_guess: "BaseFluid | None" = None,
    ):
        if mdot <= 0.0:
            raise ValueError(f"Recycle {name!r}: mdot must be > 0, got {mdot}")
        self.name = name
        self.mdot = mdot
        self._fluid_guess = fluid_guess
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def residuals(self, state: "NetworkState") -> list[float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        return [P_out - P_in, h_out - h_in]

    def fixed_node_mdot(self) -> dict[str, float]:
        return {self._outlet_node: self.mdot}

    def fluid_seed(self) -> dict[str, "BaseFluid"]:
        if self._fluid_guess is None:
            return {}
        return {self._outlet_node: self._fluid_guess}

    def recycle_fluid_delta(self, state: "NetworkState") -> float | None:
        """Max absolute per-species mass-fraction difference between the
        current guess (seeded at this component's outlet) and what actually
        propagated back around the loop to this component's inlet -- the
        tear-stream convergence check Network.solve()'s outer loop uses.

        None means "nothing to check": either fluid_guess isn't set (this
        Recycle isn't tearing a composition cycle at all), or the resolved
        fluids don't expose mass_fractions() (same duck-typed contract
        Junction._all_blendable() checks -- nothing comparable, so treated
        as trivially converged rather than guessed at, same precedent
        Junction.merge_fluids() itself follows).
        """
        if self._fluid_guess is None:
            return None
        actual = state.fluid_at(self._inlet_node)
        guess = self._fluid_guess
        if not (hasattr(guess, "mass_fractions") and hasattr(actual, "mass_fractions")):
            return None
        guess_Y = guess.mass_fractions()
        actual_Y = actual.mass_fractions()
        species = set(guess_Y) | set(actual_Y)
        if not species:
            return None
        return max(abs(actual_Y.get(s, 0.0) - guess_Y.get(s, 0.0)) for s in species)

    def update_fluid_guess(self, state: "NetworkState") -> None:
        """Replace the guess with what actually arrived at this component's
        inlet this outer iteration -- Network.solve()'s direct-substitution
        step for the next pass. No-op if fluid_guess was never set."""
        if self._fluid_guess is None:
            return
        self._fluid_guess = state.fluid_at(self._inlet_node)
