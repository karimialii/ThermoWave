from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid


class BaseComponent(ABC):
    """Interface every network component implements."""

    name: str

    @abstractmethod
    def ports(self) -> dict[str, str]:
        """Named ports ('inlet'/'outlet'/etc.) -> this component's own port id.

        Port ids are local to the component (typically derived from its name);
        Network.connect() merges two components' port ids into one shared node.
        """

    def internal_nodes(self) -> list[str]:
        """Node names this component creates internally (not shared with neighbors)."""
        return []

    @abstractmethod
    def residuals(self, state: "NetworkState") -> list[float]:
        """Residual equations contributed by this component given current state."""

    def fixed_node_values(self, fluid: "BaseFluid") -> dict[str, tuple[float, float]]:
        """node_name -> (P, h) this component fixes as a boundary condition."""
        return {}

    def fixed_node_mdot(self) -> dict[str, float]:
        """node_name -> mass flow rate [kg/s] this component fixes as a boundary condition."""
        return {}

    def guess_node_mdot(self) -> dict[str, float]:
        """node_name -> initial-guess mass flow rate [kg/s], for nodes this
        component leaves free (not in fixed_node_mdot()) rather than fixes.
        Only meaningful when a Source's own mdot is left unfixed so the
        network's total mass flow becomes a Newton unknown (e.g. closed by
        some other boundary condition, like a Sink pinning its inlet
        pressure) — without this, the solver's flat 1.0 kg/s fallback guess
        can land a map-based component's corrected mass flow so far outside
        its table that the first Jacobian is singular before Newton takes a
        step. Default: no guess (falls back to the solver's flat default).
        """
        return {}

    def free_parameters(self) -> dict[str, float]:
        """param_name -> initial guess, for scalar unknowns owned by this component
        (e.g. a compressor's shaft speed N, solved for instead of fixed, when some
        other target like PR or power is given in its place). The solver adds one
        unknown per entry and expects residuals() to contribute exactly one more
        equation per entry (typically "computed_value - target = 0"). Read back via
        state.param(f"{self.name}.<param_name>") inside residuals().

        Only used to discover which parameters are free and as a fallback
        initial guess when guess_free_parameters() below can't do better
        (e.g. no context to work from). Keep this cheap/generic; put any
        context-aware guessing in guess_free_parameters() instead.
        """
        return {}

    def guess_free_parameters(
        self, fluid: "BaseFluid", P_in: float, h_in: float, mdot: float
    ) -> dict[str, float]:
        """Like free_parameters(), but allowed to use the solver's warm-start
        guess of this component's own inlet state (and the fluid model) for a
        better initial guess than a hardcoded reference — e.g. a turbine deep
        downstream of a heater can guess its shaft speed from the actual
        (warm-started) inlet temperature it'll see there, rather than
        assuming some fixed reference temperature that's right in isolation
        but wrong wherever the component actually ends up sitting in a
        network. Must return the same keys as free_parameters(); default
        ignores the context and falls back to it.
        """
        return self.free_parameters()

    def differential_parameters(self) -> dict[str, float]:
        """param_name -> initial value, for scalar state this component
        integrates over time in Network.solve_transient() (e.g. a dynamic
        Shaft's rotor speed N, driven by inertia and net torque). The solver
        adds one unknown per entry, same as free_parameters() — but unlike
        free_parameters(), it closes it automatically rather than expecting
        a matching residual from elsewhere:
          - Network.solve() (steady state, dt=None): state_derivative()
            for that key must equal zero.
          - Network.solve_transient() (dt given): backward-Euler,
            (value - value_at_previous_step) / dt == state_derivative().
        Read back via state.param(f"{self.name}.<param_name>") inside
        residuals() and state_derivative(), same as free_parameters(). Do
        not also contribute a residual pinning this value yourself —
        the solver already supplies exactly one.
        """
        return {}

    def state_derivative(self, state: "NetworkState") -> dict[str, float]:
        """d(value)/dt for every key declared in differential_parameters(),
        evaluated at the current state. Must return exactly those keys.
        Only called for components that declare at least one differential
        parameter; the default here is never reached by such components
        (they must override it) and is only a placeholder for the rest.
        """
        return {}

    def warm_start_pairs(self) -> list[tuple[str, str]]:
        """(inlet_port, outlet_port) name pairs this component's guess_outlet()
        (or guess_outlet_for_pair(), for more than one pair) forward-
        propagates a warm-start guess across — see Solver.solve()'s guess
        propagation loop. Default: the single ("in", "out") pair, if
        this component actually has ports named that (most do); otherwise
        none, so it's silently skipped rather than erroring. Override for a
        component with differently-named or multiple flow-through port pairs
        (e.g. SimpleHeatExchanger's hot_in/hot_out and cold_in/cold_out) so
        its guess actually gets used for each of them, instead of leaving
        those downstream nodes to fall back to the solver's flat default
        guess — which, left unfixed for long enough (e.g. a map-based
        component two or three components further downstream), can put that
        component's corrected operating point so far outside its valid range
        that the very first Jacobian is singular before Newton takes a
        single step.
        """
        ports = self.ports()
        if "in" in ports and "out" in ports:
            return [("in", "out")]
        return []

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        """Rough forward estimate of this component's outlet (P, h) given a
        guessed inlet state — used only to warm-start the Newton solver's
        initial guess for downstream free nodes (see Solver._seed_x0()), not
        part of the physics itself and never appears in residuals(). Default:
        pass through unchanged, which is fine for components with a small
        P/h swing. Override it for components with a large, predictable
        swing (a Pipe's heat_loss, a Compressor/Turbine's typical PR) —
        without that, every free node in a multi-component chain starts
        Newton at whatever the network's first fixed (P, h) happens to be
        (usually an ambient Source), which can put downstream map-based
        components (e.g. a turbine well after a compressor + heater) so far
        outside their map's valid range at iteration 0 that the very first
        Jacobian is singular before Newton takes a single step.

        Only called for the single-pair case (warm_start_pairs() returning
        [("in", "out")]); a component with multiple or differently-
        named pairs implements guess_outlet_for_pair() instead — kept as a
        separate method rather than adding a pair argument here so every
        existing single-pair override keeps its original 3-argument form.
        """
        return P_in, h_in

    def guess_outlet_for_pair(
        self, pair: tuple[str, str], P_in: float, h_in: float, mdot: float
    ) -> tuple[float, float]:
        """Like guess_outlet(), for one of the pairs in a warm_start_pairs()
        list longer than the single default ("in", "out") — pair
        identifies which one is currently being propagated. Default:
        delegate to guess_outlet() (fine for the common single-pair case;
        override this instead when different pairs need different
        treatment, as SimpleHeatExchanger's hot/cold sides do).
        """
        return self.guess_outlet(P_in, h_in, mdot)

    def guess_outlet_mdot(self, pair: tuple[str, str], mdot_in: float) -> float:
        """Forward warm-start guess for a warm_start_pairs() pair's outlet
        mass flow, given a guessed inlet mdot — same role as guess_outlet()/
        guess_outlet_for_pair() but for mass flow, propagated across the
        network alongside (P, h) whenever a Source leaves mdot unfixed
        (mdot=None) and total flow is a Newton unknown instead of a given.
        Default: pass through unchanged, correct for anything that
        conserves mass across the pair (the overwhelming majority of
        components). Override for one that doesn't (e.g. SimpleCombustor
        adding fuel mass) — without that, every mass-flow unknown past such
        a component defaults to the exact same flat guess as everything
        upstream of it, which can leave two streams that are supposed to
        differ (e.g. a recuperator's hot vs. cold side once one of them has
        picked up fuel mass) landing on an artificially exact tie. That's
        harmless by itself, but if a component's own physics has a
        derivative that's only non-smooth (or only ill-conditioned) exactly
        at such a tie — e.g. an effectiveness-NTU heat exchanger's
        Cmin = min(C_hot, C_cold) — it can make the very first Jacobian
        singular before Newton takes a single step.
        """
        return mdot_in

    def fluid_flow_pairs(self) -> list[tuple[str, str]]:
        """(inlet_port, outlet_port) pairs Network._resolve_node_fluid()'s
        propagation loop uses to call outlet_fluid() -- separate from
        warm_start_pairs() (which drives the solver's unrelated (P, h)
        numeric guess) so a component can opt out of the generic single-
        inlet-per-pair fluid propagation without disturbing its own P/h
        warm-start behavior. Default: delegate to warm_start_pairs(), which
        covers every ordinary flow-through component (one inlet feeds one
        outlet, composition just passes through or is transformed by
        outlet_fluid()) with no extra work. Override to return [] for a
        genuine multi-inlet merge point that needs every inlet's resolved
        fluid at once rather than one at a time — see merge_fluids().
        """
        return self.warm_start_pairs()

    def merge_fluids(
        self, state: "NetworkState", inlet_fluids: dict[str, "BaseFluid"]
    ) -> dict[str, "BaseFluid"] | None:
        """outlet_port_name -> BaseFluid for every outlet this component can
        currently resolve, given inlet_fluids (inlet_port_name -> BaseFluid
        for every one of this component's OWN inlet ports already resolved
        this propagation pass — ports not yet resolved are simply absent).
        Default: None, meaning "not a merge point" (the overwhelming
        majority of components — anything with at most one inlet feeding
        one outlet should use outlet_fluid()/fluid_flow_pairs() instead).

        Return None if not enough inlets are resolved yet to compute a
        correct answer (e.g. only some of a Junction's several inlets are
        known this pass) — Network._resolve_node_fluid()'s fixed-point loop
        simply retries later passes, the same "not ready yet, skip and
        retry" pattern the pair-based path already uses. Only return a real
        dict once every inlet this component actually needs is present in
        inlet_fluids, since (unlike outlet_fluid()) there's no separate
        per-pair retry — get it right in one call.
        """
        return None

    def outlet_fluid(
        self, state: "NetworkState", pair: tuple[str, str], inlet_fluid: "BaseFluid"
    ) -> "BaseFluid | None":
        """BaseFluid flowing out of this component's outlet port for the
        given warm_start_pairs() pair, given the already-resolved fluid at
        its inlet port for that same pair. Default None means pass-through
        (the outlet carries the same fluid/composition as the inlet) --
        correct for the overwhelming majority of components, which don't
        change what's chemically flowing through them. Override only for a
        component whose physics genuinely changes composition (e.g.
        Combustor mixing in fuel and reacting it) to return a different
        BaseFluid instance instead.

        Called once per Newton residual evaluation (see
        Network._resolve_node_fluid()) against the network's actual current
        state -- not just a warm-start guess -- so an outlet composition
        that itself depends on a free unknown (e.g. Combustor's mdot_fuel)
        stays consistent with whatever value the solver is currently
        trying, iteration to iteration. Keep this cheap enough to call every
        residual evaluation: reuse a pre-loaded mechanism/Solution object
        rather than re-parsing one here (see Combustor's own docstring for
        why that mattered even before this hook existed), and use
        `state._cache` to memoize anything else this component's own
        residuals() would otherwise have to recompute a second time.
        """
        return None

    def report_metrics(self, state: "NetworkState") -> dict[str, float] | None:
        """Performance metrics (power, efficiency, pressure ratio, ...) for the
        components table in SolveResult.print_report(), or None (the default,
        used by Source/Sink/Pipe) to be left out of that table entirely.
        """
        return None

    def report_category(self) -> str | None:
        """Which table in SolveResult.print_report() this component's
        report_metrics() belongs in: "turbomachinery", "heat_exchanger",
        "controller", "sensor", "generator", or None (the default) to be
        left out of the categorized tables entirely (e.g. Source/Sink/Pipe,
        which also return None from report_metrics() anyway). Components
        sharing a category share one column layout, so put report_metrics()
        keys that make sense side-by-side in the same category.
        """
        return None
