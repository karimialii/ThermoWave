from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class SimpleGenerator(BaseComponent):
    """Analytic generator: fixed mechanical-to-electrical efficiency.

    Has no flow ports of its own — like Setpoint/Sensor, it's a passive
    reader, not part of the flow network. It reads a shaft component's own
    report_metrics()["power [W]"] (e.g. a Turbine's shaft power) and scales
    it down by efficiency to get electrical output:
        power_elec = shaft_power * efficiency
    Contributes zero residuals: it doesn't feed back into or constrain the
    thermodynamic solve, only reports a derived quantity from it. A map-
    based Generator (speed-vs-torque characteristic) is available as a
    separate, more detailed component.
    """

    def __init__(self, name: str, component: BaseComponent, efficiency: float):
        self.name = name
        self.component = component
        self.efficiency = efficiency

    def ports(self) -> dict[str, str]:
        return {}

    def report_category(self) -> str:
        return "generator"

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        shaft_metrics = self.component.report_metrics(state)
        if shaft_metrics is None or "power [W]" not in shaft_metrics:
            raise ValueError(
                f"SimpleGenerator {self.name!r} reads shaft power from "
                f"{self.component.name!r}, but its report_metrics() doesn't "
                f"expose 'power [W]' (got: {sorted(shaft_metrics) if shaft_metrics else []})"
            )
        shaft_power = shaft_metrics["power [W]"]
        metrics = {"power [W]": shaft_power * self.efficiency, "eta [-]": self.efficiency}
        if "N [rev/min]" in shaft_metrics:
            metrics["N [rev/min]"] = shaft_metrics["N [rev/min]"]
        return metrics
