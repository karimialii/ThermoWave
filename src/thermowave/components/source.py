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
    """

    def __init__(
        self, name: str, P: float, T: float, mdot: float | None, mdot_guess: float = 1.0
    ):
        self.name = name
        self.P_si = settings.pressure_to_si(P)
        self.T_si = settings.temperature_to_si(T)
        self.mdot = mdot
        self.mdot_guess = mdot_guess
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"out": self._outlet_node}

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def fixed_node_values(self, fluid: "BaseFluid") -> dict[str, tuple[float, float]]:
        h = fluid.enthalpy_pt(self.P_si, self.T_si)
        return {self._outlet_node: (self.P_si, h)}

    def fixed_node_mdot(self) -> dict[str, float]:
        if self.mdot is None:
            return {}
        return {self._outlet_node: self.mdot}

    def guess_node_mdot(self) -> dict[str, float]:
        if self.mdot is None:
            return {self._outlet_node: self.mdot_guess}
        return {}
