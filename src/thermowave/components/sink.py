from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class Sink(BaseComponent):
    """Open boundary terminating a network branch.

    Leave P=None (default) to fix no state at all — a pass-through exit that
    just absorbs whatever mass flow already gets fixed elsewhere (typically
    a Source's own mdot). Give P to instead pin the inlet node's pressure to
    that fixed value, via one residual (P_in - P), leaving h free (set by
    whatever upstream energy balance reaches this node): this is the
    physical closure a real open exhaust actually has (it expands back down
    to ambient), and pairs with Source(mdot=None) — total mass flow is then
    a Newton unknown, solved for as whatever value is self-consistent with
    the shaft speed and every component's own characteristic in between,
    rather than an externally dictated input. Adds exactly the one equation
    that freeing mdot removes, so the system stays square.
    """

    def __init__(self, name: str, P: float | None = None):
        self.name = name
        self._inlet_node = f"{name}.in"
        if P is None:
            self.P_si = None
        else:
            from thermowave.core.settings import settings

            self.P_si = settings.pressure_to_si(P)

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node}

    def residuals(self, state: "NetworkState") -> list[float]:
        if self.P_si is None:
            return []
        P_in, _ = state.node(self._inlet_node)
        return [P_in - self.P_si]
