from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.components.sensor import SENSOR_QUANTITIES as _ALLOWED_QUANTITIES

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class Controller(BaseComponent):
    """Drives a component's free parameter until a Sensor reads a target value.

    Has no ports of its own — it contributes exactly one residual:
    sensor.report_metrics(state)[quantity] - value == 0. Same closed-loop
    idea as Setpoint (a target instead of a direct input, pinned down by
    leaning on a free parameter the target component already declared), but
    where Setpoint reads the *target component's own* report_metrics()
    (e.g. a compressor's own power/PR/eta_s), Controller reads an
    independent Sensor's measurement instead — e.g. drive a Compressor's
    free N until a Sensor sitting on a downstream pipe outlet reads a
    target temperature, mirroring a real control loop where the measured
    and actuated quantities live in different places in the plant.

    Raises ValueError at construction if the target component doesn't
    currently declare free_param as free, for the same reason as Setpoint:
    failing fast here beats a mismatched Newton system discovered later.
    """

    def __init__(
        self,
        name: str,
        sensor: BaseComponent,
        quantity: str,
        component: BaseComponent,
        free_param: str,
        value: float,
    ):
        if free_param not in component.free_parameters():
            raise ValueError(
                f"Controller {name!r} actuates {component.name!r}.{free_param}, but "
                f"that component doesn't currently declare {free_param!r} as free — "
                f"pass None for it at construction so it becomes a solvable unknown."
            )

        self.name = name
        self.sensor = sensor
        self.quantity = quantity
        self.component = component
        self.free_param = free_param
        self.value = value

    def ports(self) -> dict[str, str]:
        return {}

    def report_category(self) -> str:
        return "controller"

    def residuals(self, state: "NetworkState") -> list[float]:
        metrics = self.sensor.report_metrics(state)
        if self.quantity not in metrics:
            raise ValueError(
                f"Controller {self.name!r} reads quantity {self.quantity!r} from "
                f"sensor {self.sensor.name!r}, but it doesn't expose that reading "
                f"(got: {sorted(metrics)}; try one of {_ALLOWED_QUANTITIES})"
            )
        return [metrics[self.quantity] - self.value]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        measured = self.sensor.report_metrics(state)[self.quantity]
        return {
            "target [-]": self.value,
            "measured [-]": measured,
            "error [-]": measured - self.value,
        }
