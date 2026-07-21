from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

# The quantities Sensor exposes via report_metrics() — also the set
# Controller/PIDController accept for their own `quantity` arg, since both
# read a Sensor (or Sensor-shaped) component the same way. Defined once
# here rather than duplicated in each of those files.
SENSOR_QUANTITIES = ("P [Pa]", "T [K]", "h [J/kg]", "mdot [kg/s]")


class Sensor(BaseComponent):
    """Passive measurement tap on a network node.

    Has a single port ("tap") — connect() it to whatever node you want to
    read (e.g. connect(pipe, "out", sensor, "tap")), which merges the
    sensor into that existing node rather than creating a new one. It
    contributes zero residuals and fixes nothing: purely a read, never a
    perturbation, so wiring one in anywhere never changes the solved state
    or adds an unknown.

    report_metrics() exposes the reading — "P [Pa]", "T [K]", "h [J/kg]",
    and "mdot [kg/s]" (omitted if this node's mdot isn't part of the solve,
    e.g. a node with no flow through it). A Controller reads these the same
    way Setpoint reads a component's own report_metrics().
    """

    def __init__(self, name: str):
        self.name = name
        self._node = f"{name}.tap"

    def ports(self) -> dict[str, str]:
        return {"tap": self._node}

    def report_category(self) -> str:
        return "sensor"

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P, h = state.node(self._node)
        metrics = {"P [Pa]": P, "h [J/kg]": h, "T [K]": state.fluid_at(self._node).temperature_ph(P, h)}
        try:
            metrics["mdot [kg/s]"] = state.mdot(self._node)
        except KeyError:
            pass
        return metrics
