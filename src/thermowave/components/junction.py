from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.fluids.cantera_fluid import _CanteraCompositionFluid

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid

_MDOT_FLOOR = 1.0e-9  # kg/s, guards h_mix's division against a transient
# all-zero-inflow Newton iterate — see residuals()'s docstring note below.


class Junction(BaseComponent):
    """Splitter/Merger: any number of inlets combine and split into any
    number of outlets.

    Modeled as a zero-loss, zero-volume plenum: every outlet shares one
    common pressure (no pressure drop across the junction) and carries the
    same mass-weighted mixture enthalpy of all the inlets (perfect,
    instantaneous mixing). Outlet mass flow is a fixed fraction of the total
    inlet mass flow — split_fractions is a modeling input, not something
    solved from downstream flow resistance, so it must be supplied by the
    caller and should sum to 1.

    Junction only produces residuals for its own outlets: inlet nodes are
    owned (and already fully constrained) by whichever upstream component
    feeds them, exactly like every other component's inlet port. The
    junction's common reference pressure is taken from the first inlet;
    with multiple inlets at genuinely different pressures this is a
    simplification (no attempt is made to solve for a physically
    self-consistent merge pressure).

    n_inlets == 1 with n_outlets > 1 is a pure Splitter; n_inlets > 1 with
    n_outlets == 1 is a pure Merger (split_fractions is then just [1.0]).

    Composition mixing (merge_fluids()): if every inlet's fluid (see
    NetworkState.fluid_at()) is already the same object — the overwhelming
    majority of networks, where nothing upstream changes composition — that
    fluid just passes straight through to every outlet, no blending needed.
    If inlets genuinely differ *and* every one of them exposes a
    mass_fractions() method and a `mechanism` attribute on a common
    mechanism (both CanteraFluid and Combustor's own equilibrium-product
    fluid do — duck-typed rather than an isinstance() check against one
    specific class, precisely so two *different* combustor exhaust streams
    merging, the realistic case this exists for, actually qualifies), this
    computes a real mass-weighted mixture (each inlet's own fixed mass
    fractions, weighted by its mdot, summed and renormalized) and exposes
    it as a fresh fluid view for every outlet — a genuine merge, not just a
    label. Anything else (inlets missing that contract, or on different
    mechanisms — no common species set to blend across) falls back to the
    first inlet's fluid: an explicit, documented simplification rather than
    a real merge, the same class of scope boundary Combustor's own
    outlet_fluid() gating makes for the same underlying reason (mixing
    fluid models that don't share Cantera's absolute enthalpy datum would
    silently corrupt any downstream energy balance).
    """

    def __init__(
        self,
        name: str,
        n_inlets: int,
        n_outlets: int,
        split_fractions: list[float] | None = None,
    ):
        if split_fractions is None:
            split_fractions = [1.0 / n_outlets] * n_outlets
        if len(split_fractions) != n_outlets:
            raise ValueError(
                f"split_fractions has {len(split_fractions)} entries, "
                f"expected n_outlets={n_outlets}"
            )

        self.name = name
        self.n_inlets = n_inlets
        self.n_outlets = n_outlets
        self.split_fractions = split_fractions
        self._inlet_nodes = [f"{name}.in{i}" for i in range(n_inlets)]
        self._outlet_nodes = [f"{name}.out{i}" for i in range(n_outlets)]
        # Lazily created the first time merge_fluids() actually needs to
        # blend two differently-composed CanteraFluid inlets -- most
        # Junctions never touch Cantera at all, so importing/loading it
        # eagerly here would be pure overhead for those.
        self._mix_gas = None

    def ports(self) -> dict[str, str]:
        ports = {f"in{i}": node for i, node in enumerate(self._inlet_nodes)}
        ports.update({f"out{i}": node for i, node in enumerate(self._outlet_nodes)})
        return ports

    def fluid_flow_pairs(self) -> list[tuple[str, str]]:
        # Junction resolves every outlet's fluid entirely through
        # merge_fluids() (it needs every inlet at once, not one pair at a
        # time) -- opting out of the generic pair-based path here stops it
        # from racing ahead and locking in a wrong "only in0 matters"
        # answer before merge_fluids() has every inlet it actually needs.
        return []

    def merge_fluids(
        self, state: "NetworkState", inlet_fluids: dict[str, "BaseFluid"]
    ) -> dict[str, "BaseFluid"] | None:
        resolved = [inlet_fluids.get(f"in{i}") for i in range(self.n_inlets)]
        if any(f is None for f in resolved):
            return None  # not every inlet known yet this propagation pass

        first = resolved[0]
        if all(f is first for f in resolved):
            mixed = first
        elif self._all_blendable(resolved):
            mixed = self._blend(state, resolved)
        else:
            # Can't do a principled mass-weighted composition blend across
            # fluid models that don't all share a common composition/
            # mechanism -- see this class's own docstring. Falls back to
            # the first inlet's fluid rather than raising or guessing.
            mixed = first

        return {f"out{i}": mixed for i in range(self.n_outlets)}

    @staticmethod
    def _all_blendable(fluids: list["BaseFluid"]) -> bool:
        first = fluids[0]
        if not (hasattr(first, "mass_fractions") and hasattr(first, "mechanism")):
            return False
        return all(
            hasattr(f, "mass_fractions")
            and hasattr(f, "mechanism")
            and f.mechanism == first.mechanism
            for f in fluids
        )

    def _blend(self, state: "NetworkState", fluids: list["BaseFluid"]) -> "BaseFluid":
        """Mass-weighted mixture of two or more same-mechanism inlets whose
        fluid exposes mass_fractions() (see this class's own docstring for
        why that's a duck-typed check, not isinstance()) -- each inlet's own
        fixed mass fractions, weighted by its mdot, summed and renormalized.
        Reuses one lazily-created Solution (self._mix_gas) across every call
        instead of re-parsing the mechanism per Newton residual evaluation,
        the same cost-avoidance Combustor's own _product_gas already
        established.
        """
        mdots = [state.mdot(node) for node in self._inlet_nodes]
        mdot_total = sum(mdots)
        if mdot_total < _MDOT_FLOOR:
            return fluids[0]  # degenerate all-zero-flow guess; any answer is as good as any

        species_mass: dict[str, float] = {}
        for fluid, mdot in zip(fluids, mdots):
            for species, y in fluid.mass_fractions().items():
                species_mass[species] = species_mass.get(species, 0.0) + mdot * y
        blended_Y = {species: mass / mdot_total for species, mass in species_mass.items()}

        if self._mix_gas is None:
            import cantera as ct

            self._mix_gas = ct.Solution(fluids[0].mechanism)
        return _CanteraCompositionFluid(f"{self.name}.mix", self._mix_gas, blended_Y, fluids[0].mechanism)

    def warm_start_pairs(self) -> list[tuple[str, str]]:
        # Port names here are in0/out0/... rather than the single-pair
        # default's in/out, so BaseComponent's default warm_start_pairs()
        # (which only looks for a port literally named "in"/"out") silently
        # skips this component — every downstream node then falls back to
        # the solver's flat default guess instead of a propagated one. Pair
        # the first inlet (the same one residuals() below already treats as
        # the reference for P_ref) with every outlet; good enough for a
        # zero-loss, zero-drop plenum where every outlet's guess only needs
        # to be in the right order of magnitude, not exact.
        if not self._inlet_nodes or not self._outlet_nodes:
            return []
        return [("in0", f"out{i}") for i in range(self.n_outlets)]

    def residuals(self, state: "NetworkState") -> list[float]:
        inlet_mdots = [state.mdot(node) for node in self._inlet_nodes]
        inlet_states = [state.node(node) for node in self._inlet_nodes]
        mdot_total = sum(inlet_mdots)
        # A genuinely all-zero-inflow state (a plausible mid-iteration
        # Newton guess, or a transient's t=0 before flow ramps up) would
        # otherwise raise ZeroDivisionError here instead of letting the
        # solver's own step-limiting/damping machinery recover — h_mix is
        # physically undefined at mdot_total==0 either way, so flooring the
        # denominator just makes it 0/floor == 0 instead of NaN/crash.
        mdot_total_safe = mdot_total if abs(mdot_total) > _MDOT_FLOOR else _MDOT_FLOOR
        h_mix = (
            sum(mdot * h for mdot, (_, h) in zip(inlet_mdots, inlet_states)) / mdot_total_safe
        )
        P_ref = inlet_states[0][0]

        out: list[float] = []
        for i, node in enumerate(self._outlet_nodes):
            P_out, h_out = state.node(node)
            out.append(P_out - P_ref)
            out.append(h_out - h_mix)
            out.append(state.mdot(node) - self.split_fractions[i] * mdot_total)
        return out
