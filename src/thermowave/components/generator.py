from __future__ import annotations

import math
from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.maps.torque_speed_map import TorqueSpeedMap

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_RAD_PER_MIN_TO_RAD_PER_S = 2.0 * math.pi / 60.0


class Generator(BaseComponent):
    """Generator driven by a speed-vs-torque characteristic map.

    Has no flow ports of its own — like SimpleGenerator, it's a passive
    reader of a shaft component's own report_metrics()["N [rev/min]"] (e.g.
    a Turbine's shaft speed). Mechanical power available at that speed comes
    from the map's own torque curve rather than from the shaft component's
    reported power:
        omega = N * 2*pi/60
        power_mech = map.torque(N) * omega
        power_elec = power_mech * efficiency
    (efficiency defaults to 1.0, i.e. the map's torque curve is assumed to
    already be the generator's electrical output rating; pass a lower value
    if the map is a mechanical/shaft rating and electrical losses still need
    to be applied on top of it.) Contributes zero residuals: it doesn't feed
    back into or constrain the thermodynamic solve, only reports a derived
    quantity from it.
    """

    def __init__(
        self,
        name: str,
        component: BaseComponent,
        map_path: str,
        efficiency: float = 1.0,
    ):
        self.name = name
        self.component = component
        self.map = TorqueSpeedMap.from_file(map_path)
        self.efficiency = efficiency

    def ports(self) -> dict[str, str]:
        return {}

    def report_category(self) -> str:
        return "generator"

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        shaft_metrics = self.component.report_metrics(state)
        if shaft_metrics is None or "N [rev/min]" not in shaft_metrics:
            raise ValueError(
                f"Generator {self.name!r} reads shaft speed from "
                f"{self.component.name!r}, but its report_metrics() doesn't "
                f"expose 'N [rev/min]' (got: {sorted(shaft_metrics) if shaft_metrics else []})"
            )
        N = shaft_metrics["N [rev/min]"]
        omega = N * _RAD_PER_MIN_TO_RAD_PER_S
        torque = self.map.torque(N)
        power_mech = torque * omega
        return {
            "power [W]": power_mech * self.efficiency,
            "eta [-]": self.efficiency,
            "N [rev/min]": N,
        }
