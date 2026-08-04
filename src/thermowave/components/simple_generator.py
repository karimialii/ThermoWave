from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class SimpleGenerator(BaseComponent):
    """Analytic generator: fixed mechanical-to-electrical efficiency.

    Has no flow ports of its own — like Setpoint/Sensor, it's a passive
    reader, not part of the flow network. It reads shaft power through a
    "power" signal port (e.g. a Turbine's shaft power), wired with
    network.connect(gen, "power", turbine, "power", kind="signal"), and
    scales it down by efficiency to get electrical output:
        power_elec = shaft_power * efficiency
    Contributes zero residuals: it doesn't feed back into or constrain the
    thermodynamic solve, only reports a derived quantity from it. It also
    exposes its own "shaft" mechanical port -- connect it to the same
    machine's "shaft" port (kind="mechanical") so N shows up in this
    component's own report_metrics() too; leave it unconnected if that's
    not needed. A map-based Generator (speed-vs-torque characteristic) is
    available as a separate, more detailed component.
    """

    def __init__(self, name: str, efficiency: float):
        self.name = name
        self.efficiency = efficiency
        self._shaft_port = f"{name}.shaft"
        self._power_port = f"{name}.power"

    def ports(self) -> dict[str, str]:
        return {}

    def mechanical_ports(self) -> dict[str, str]:
        return {"shaft": self._shaft_port}

    def signal_ports(self) -> dict[str, str]:
        return {"power": self._power_port}

    def report_category(self) -> str:
        return "generator"

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        shaft_power = state.signal(self._power_port)
        metrics = {"power [W]": shaft_power * self.efficiency, "eta [-]": self.efficiency}
        try:
            metrics["N [rev/min]"] = state.N(self._shaft_port)
        except KeyError:
            pass
        return metrics
