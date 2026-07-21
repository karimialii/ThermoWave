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
    required report_metrics()["power [W]"] (e.g. an electrically-driven
    Compressor's own shaft power demand, with no Turbine on its shaft to
    supply it) and reports the electrical power that must be drawn to
    supply it, given motor efficiency:
        power_elec = shaft_power_required / efficiency
    the reverse of SimpleGenerator's power_elec = shaft_power * efficiency.

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

    def __init__(self, name: str, component: BaseComponent, efficiency: float):
        if not (0.0 < efficiency <= 1.0):
            raise ValueError(
                f"ElectricMotor {name!r}: efficiency must be in (0, 1], got {efficiency}"
            )
        self.name = name
        self.component = component
        self.efficiency = efficiency

    def ports(self) -> dict[str, str]:
        return {}

    def report_category(self) -> str:
        return "motor"

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        mech_metrics = self.component.report_metrics(state)
        if mech_metrics is None or "power [W]" not in mech_metrics:
            raise ValueError(
                f"ElectricMotor {self.name!r} reads required shaft power from "
                f"{self.component.name!r}, but its report_metrics() doesn't "
                f"expose 'power [W]' (got: {sorted(mech_metrics) if mech_metrics else []})"
            )
        shaft_power = mech_metrics["power [W]"]
        metrics = {"power [W]": shaft_power / self.efficiency, "eta [-]": self.efficiency}
        if "N [rev/min]" in mech_metrics:
            metrics["N [rev/min]"] = mech_metrics["N [rev/min]"]
        return metrics
