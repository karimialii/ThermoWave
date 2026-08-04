from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.components.sensor import SENSOR_QUANTITIES as _ALLOWED_QUANTITIES
from thermowave.core.network import TargetValue

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class Controller(BaseComponent):
    """Drives a component's free parameter until a Sensor reads a target value.

    Has no ports of its own — it contributes exactly one residual:
    sensor.report_metrics(state)[quantity] - value(state) == 0. Same
    closed-loop idea as Setpoint (a target instead of a direct input,
    pinned down by leaning on a free parameter the target component
    already declared), but where Setpoint reads the *target component's
    own* report_metrics() (e.g. a compressor's own power/PR/eta_s),
    Controller reads an independent Sensor's measurement instead — e.g.
    drive a Compressor's free N until a Sensor sitting on a downstream
    pipe outlet reads a target temperature, mirroring a real control loop
    where the measured and actuated quantities live in different places in
    the plant.

    value: a plain float (the common case), or a callable state -> float
    re-evaluated fresh every residual call — the same widening Setpoint
    has, and for the same reason: tying this Sensor's reading to another
    live quantity (another Sensor, or any component's report_metrics())
    instead of a constant, e.g. matching two streams' temperatures at a
    merge point:
        Controller(
            name="match_T", sensor=sensor_a, quantity="T [K]",
            component=recompressor, free_param="N",
            value=lambda s: sensor_b.report_metrics(s)["T [K]"],
        )
    See Setpoint's docstring for the fuller rationale — both classes share
    this behavior via the same TargetValue type.

    free_param may also name a mechanical port left free (e.g. a Turbine/
    Compressor with N=None exposes "shaft") instead of an ordinary free
    parameter -- see Setpoint's docstring for the same widening.

    Raises ValueError at construction if the target component doesn't
    currently declare free_param as free (in either sense), for the same
    reason as Setpoint: failing fast here beats a mismatched Newton system
    discovered later.
    """

    def __init__(
        self,
        name: str,
        sensor: BaseComponent,
        quantity: str,
        component: BaseComponent,
        free_param: str,
        value: TargetValue,
    ):
        mech_ports = component.mechanical_ports()
        self._mechanical = (
            free_param in mech_ports and mech_ports[free_param] in component.free_mechanical_ports()
        )
        if not self._mechanical and free_param not in component.free_parameters():
            raise ValueError(
                f"Controller {name!r} actuates {component.name!r}.{free_param}, but "
                f"that component doesn't currently declare {free_param!r} as a free "
                f"parameter or a free mechanical port — pass None for it at "
                f"construction so it becomes a solvable unknown."
            )

        self.name = name
        self.sensor = sensor
        self.quantity = quantity
        self.component = component
        self.free_param = free_param
        self.value = value

    def _target(self, state: "NetworkState") -> float:
        return self.value(state) if callable(self.value) else self.value

    def ports(self) -> dict[str, str]:
        return {}

    def closes_parameters(self) -> list[str]:
        if self._mechanical:
            return []
        return [f"{self.component.name}.{self.free_param}"]

    def closes_mechanical_nodes(self) -> list[str]:
        network = getattr(self, "_network", None)
        if not self._mechanical or network is None:
            return []
        port_id = self.component.mechanical_ports()[self.free_param]
        return [network._canonical(port_id)]

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
        return [metrics[self.quantity] - self._target(state)]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        measured = self.sensor.report_metrics(state)[self.quantity]
        target = self._target(state)
        return {
            "target [-]": target,
            "measured [-]": measured,
            "error [-]": measured - target,
        }
