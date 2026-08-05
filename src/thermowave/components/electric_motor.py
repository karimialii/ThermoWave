from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class ElectricMotor(BaseComponent):
    """Electrically-driven mechanical power source — the inverse of
    SimpleGenerator/Generator.

    Has no flow ports of its own: like those, it's a passive reader, not
    part of the flow network. It reads a mechanical component's own
    required shaft power through a "power" signal port (e.g. an
    electrically-driven Compressor's own shaft power demand, with no
    Turbine on its shaft to supply it), wired with
    network.connect(motor, "power", compressor, "power", kind="signal"),
    and reports the electrical power that must be drawn to supply it, given
    motor efficiency:
        power_elec = shaft_power_required / efficiency
    the reverse of SimpleGenerator's power_elec = shaft_power * efficiency.
    It also exposes its own "shaft" mechanical port -- connect it to the
    same machine's "shaft" port (kind="mechanical") so N shows up in this
    component's own report_metrics() too; leave it unconnected if that's
    not needed.

    Contributes zero residuals: it doesn't feed back into or constrain the
    thermodynamic solve, only reports a derived quantity from it — the
    mechanical component's own free speed, if it has one (e.g. a map-based
    Compressor's N=None), still needs its own Setpoint/Controller to pin
    down, exactly as it would with no motor present at all. ElectricMotor
    only answers "how much electricity does driving this cost," the same
    narrow role SimpleGenerator plays for "how much electricity does this
    produce" — for an electrically-driven compressor/pump with genuinely no
    other shaft input, not a general two-way mechanical coupling (see Shaft
    for that).
    """

    def __init__(self, name: str, efficiency: float):
        if not (0.0 < efficiency <= 1.0):
            raise ValueError(
                f"ElectricMotor {name!r}: efficiency must be in (0, 1], got {efficiency}"
            )
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

    def shaft_sign(self) -> float:
        return 1.0  # delivers power to the shaft (electrical power in,
        # scaled down by efficiency, drives the shaft)

    def report_category(self) -> str:
        return "motor"

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        shaft_power = state.signal(self._power_port)
        metrics = {"power [W]": shaft_power / self.efficiency, "eta [-]": self.efficiency}
        try:
            metrics["N [rev/min]"] = state.N(self._shaft_port)
        except KeyError:
            pass
        return metrics
