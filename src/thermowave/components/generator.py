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

    Has no flow ports of its own -- like SimpleGenerator, it's a passive
    reader of shaft speed, but reads it through a real "shaft" mechanical
    port instead of holding a direct object reference: wire it with
    network.connect(generator, "shaft", turbine, "shaft", kind="mechanical").
    Mechanical power available at that speed comes from the map's own
    torque curve rather than from the shaft component's reported power:
        omega = N * 2*pi/60
        power_mech = map.torque(N) * omega
        power_elec = power_mech * efficiency
    (efficiency defaults to 1.0, i.e. the map's torque curve is assumed to
    already be the generator's electrical output rating; pass a lower value
    if the map is a mechanical/shaft rating and electrical losses still need
    to be applied on top of it.) Contributes zero residuals: it doesn't feed
    back into or constrain the thermodynamic solve, only reports a derived
    quantity from it. Also exposes a "power" signal port publishing the same
    electrical output, for anything downstream that wants to read it via
    kind="signal" instead of report_metrics().
    """

    def __init__(
        self,
        name: str,
        map_path: str,
        efficiency: float = 1.0,
    ):
        if not (0.0 < efficiency <= 1.0):
            raise ValueError(
                f"Generator {name!r}: efficiency must be in (0, 1], got {efficiency}"
            )
        self.name = name
        self.map = TorqueSpeedMap.from_file(map_path)
        self.efficiency = efficiency
        self._shaft_port = f"{name}.shaft"
        self._power_port = f"{name}.power"

    def ports(self) -> dict[str, str]:
        return {}

    def mechanical_ports(self) -> dict[str, str]:
        return {"shaft": self._shaft_port}

    def signal_ports(self) -> dict[str, str]:
        return {"power": self._power_port}

    def shaft_sign(self) -> float:
        return -1.0  # draws power from the shaft to generate electricity

    def report_category(self) -> str:
        return "generator"

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def _power(self, state: "NetworkState") -> tuple[float, float]:
        N = state.N(self._shaft_port)
        omega = N * _RAD_PER_MIN_TO_RAD_PER_S
        torque = self.map.torque(N)
        power_mech = torque * omega
        return N, power_mech * self.efficiency

    def provided_signal_values(self, state: "NetworkState") -> dict[str, float]:
        _, power_elec = self._power(state)
        return {self._power_port: power_elec}

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        N, power_elec = self._power(state)
        return {
            "power [W]": power_elec,
            "eta [-]": self.efficiency,
            "N [rev/min]": N,
        }
