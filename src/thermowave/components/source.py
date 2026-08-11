from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.core.settings import settings

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid


class Source(BaseComponent):
    """Boundary condition fixing outlet pressure and temperature, and
    (usually) mass flow rate.

    Leave mdot=None to leave mass flow unfixed instead — it then becomes a
    Newton unknown like any other free node mdot, for networks where total
    flow is meant to come out of the solve rather than be dictated up front
    (e.g. a turbomachinery loop closed by a Sink pinning exit pressure
    instead: the mass flow that actually reaches that pressure is
    whatever the compressor/turbine maps imply at the given shaft speed,
    not a value picked in advance). mdot_guess seeds that unknown's initial
    guess (see BaseComponent.guess_node_mdot); pick something in the right
    order of magnitude for whatever's downstream — a poor guess can put a
    map-based component so far outside its table that the first Jacobian is
    singular. Ignored (no free unknown to guess) when mdot is given.

    fluid: give a BaseFluid here to seed this Source's own outlet with a
    DIFFERENT fluid than the Network's own default (e.g. a secondary,
    non-mixing loop -- an oil/HTF stream feeding one side of a
    HeatExchanger whose other side carries the network's usual working
    fluid). Leave at the default None for the ordinary case (this Source's
    fluid is just the network's own fluid, the only behavior that existed
    before this parameter did).

    Every other component (HeatExchanger among them) already reads each of
    its own ports' fluid independently via NetworkState.fluid_at(), so a
    second fluid introduced this way flows through the rest of the network
    correctly with no further changes -- what's genuinely new here is only
    the seed itself: fixed_node_values() below uses `fluid` (if given)
    instead of the network-wide default to compute this Source's own (P, h),
    and fluid_seed() (see BaseComponent's own docstring) propagates that same
    fluid forward from this Source's outlet, the same mechanism Recycle's own
    `fluid_guess` already uses for composition cycles -- Network._resolve_
    node_fluid()'s "does anything in this network change composition" fast-
    path check also had to learn to look for fluid_seed(), not just
    outlet_fluid()/merge_fluids() (see Network._fluid_propagation_plan()),
    since a plain non-recirculating second-fluid branch (Source -> ... ->
    Sink, no Combustor/Junction anywhere) previously had no component that
    fast-path check recognized as composition-changing at all -- every node
    would have silently resolved to the network's default fluid instead.
    """

    def __init__(
        self,
        name: str,
        P: float,
        T: float,
        mdot: float | None,
        mdot_guess: float = 1.0,
        fluid: "BaseFluid | None" = None,
    ):
        if P <= 0.0:
            raise ValueError(f"Source {name!r}: P must be > 0, got {P}")
        T_si = settings.temperature_to_si(T)
        if T_si <= 0.0:
            raise ValueError(f"Source {name!r}: T must be > 0 absolute, got {T}")
        if mdot is not None and mdot < 0.0:
            raise ValueError(f"Source {name!r}: mdot must be >= 0, got {mdot}")

        self.name = name
        self.P_si = settings.pressure_to_si(P)
        self.T_si = T_si
        self.mdot = mdot
        self.mdot_guess = mdot_guess
        self.fluid = fluid
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"out": self._outlet_node}

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def fixed_node_values(self, fluid: "BaseFluid") -> dict[str, tuple[float, float]]:
        f = self.fluid if self.fluid is not None else fluid
        h = f.enthalpy_pt(self.P_si, self.T_si)
        return {self._outlet_node: (self.P_si, h)}

    def fixed_node_mdot(self) -> dict[str, float]:
        if self.mdot is None:
            return {}
        return {self._outlet_node: self.mdot}

    def guess_node_mdot(self) -> dict[str, float]:
        if self.mdot is None:
            return {self._outlet_node: self.mdot_guess}
        return {}

    def fluid_seed(self) -> dict[str, "BaseFluid"]:
        if self.fluid is None:
            return {}
        return {self._outlet_node: self.fluid}
