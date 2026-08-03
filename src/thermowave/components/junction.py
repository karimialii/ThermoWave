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
    instantaneous mixing). Outlet mass flow is a fraction of the total inlet
    mass flow — each entry of split_fractions is either a fixed number you
    supply, or None to leave it a free Newton unknown instead (the fixed
    ones plus every free one but the last must sum to 1, the same
    "otherwise unconstrained value" role Compressor's N=None plays for
    shaft speed). This is exactly right for a tree (the correct fraction is
    knowable in advance from topology and demand alone) and exactly wrong
    for a loop (where the split between two parallel paths is only known
    once both paths' resistances are, an unknown a real network solver has
    to determine) — free entries exist for the latter case, where
    something else in the network (a Setpoint tying this Junction's own
    split to a value, or a Controller/second Setpoint tying some other
    component's reading to another) supplies the closing equation(s) that
    pin them down; see free_parameters()'s docstring below for exactly how
    many and which.

    split_fractions omitted entirely (the default) still means an equal,
    fully fixed split — unchanged from before free entries existed.

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
        split_fractions: list[float | None] | None = None,
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
        # Every None entry is a free unknown except the last one, which is
        # instead computed each residual call as 1 - sum(every other entry)
        # -- see free_parameters()'s docstring for why one is deliberately
        # excluded rather than adding an (n_outlets)-th "sums to 1" residual
        # for it.
        none_indices = [i for i, f in enumerate(split_fractions) if f is None]
        self._free_indices = none_indices[:-1]
        self._dependent_index = none_indices[-1] if none_indices else None
        self._inlet_nodes = [f"{name}.in{i}" for i in range(n_inlets)]
        self._outlet_nodes = [f"{name}.out{i}" for i in range(n_outlets)]
        # Lazily created the first time merge_fluids() actually needs to
        # blend two differently-composed CanteraFluid inlets -- most
        # Junctions never touch Cantera at all, so importing/loading it
        # eagerly here would be pure overhead for those.
        self._mix_gas = None

    def free_parameters(self) -> dict[str, float]:
        """One entry per None in split_fractions, except the last such
        entry (see __init__'s docstring: that one is instead computed from
        the others each residual call, not a separate Newton unknown).
        Named f"frac{i}" for outlet index i, e.g. a 2-outlet Junction with
        split_fractions=[None, None] declares exactly one free parameter,
        "frac0" — outlet 1's fraction follows as 1 - frac0, needing no
        residual of its own. With k free entries overall, this declares
        k - 1 unknowns, each needing its own closing equation from
        elsewhere in the network (a Setpoint/Controller, typically); the
        k-th is closed automatically, by construction, not by anything
        external.
        """
        return {f"frac{i}": 1.0 / self.n_outlets for i in self._free_indices}

    def _effective_split_fractions(self, state: "NetworkState") -> list[float]:
        fracs = list(self.split_fractions)
        for i in self._free_indices:
            fracs[i] = state.param(f"{self.name}.frac{i}")
        if self._dependent_index is not None:
            fracs[self._dependent_index] = 1.0 - sum(
                f for i, f in enumerate(fracs) if i != self._dependent_index
            )
        return fracs

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
        fracs = self._effective_split_fractions(state)

        out: list[float] = []
        for i, node in enumerate(self._outlet_nodes):
            P_out, h_out = state.node(node)
            out.append(P_out - P_ref)
            out.append(h_out - h_mix)
            out.append(state.mdot(node) - fracs[i] * mdot_total)
        return out
