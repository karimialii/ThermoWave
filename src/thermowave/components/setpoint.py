from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.core.network import TargetValue

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class Setpoint(BaseComponent):
    """Drives one of a component's report_metrics() outputs to a target value
    by leaning on a free parameter that component already declared.

    Has no ports of its own — it contributes exactly one residual:
    component.report_metrics(state)[target_metric] - value(state) == 0. It
    doesn't "free" anything itself; the target component must already
    expose the named parameter as free (e.g. Compressor(..., N=None)),
    which is where the matching extra Newton unknown comes from. Setpoint
    just supplies the other half: the equation that pins that unknown down.

    This generalizes "give a target instead of a direct input" to any
    component/metric pair without teaching every component a bespoke
    N-or-PR-or-power constructor: e.g. tie a Compressor's free N to a target
    power, tie it to a target PR, or (once other components expose their own
    free_parameters()) tie a Valve's opening to a target downstream pressure.

    value: a plain float (the common case — a fixed design target), or a
    callable state -> float re-evaluated fresh every residual call, for
    tying the target metric to something that's itself part of the solve
    rather than a constant — e.g. matching two streams' temperatures at a
    merge point (a real design rule in recuperated/recompression cycles,
    where mixing streams at different temperatures destroys exergy):
        Setpoint(
            name="match_T", component=recompressor, free_param="N",
            target_metric="T_out [K]",
            value=lambda s: sensor_on_other_stream.report_metrics(s)["T [K]"],
        )
    (Sensor — thermowave.components.sensor — is the usual way to read a
    plain node's state for this; report_metrics() on any component works
    too.) A callable value has no closing residual of its own to supply —
    it's read like any other live network quantity, the same way
    target_metric already is; whatever free parameter it happens to depend
    on (if any) must be closed elsewhere, exactly as for a constant target.

    free_param may also name a mechanical port left free (e.g. a Turbine/
    Compressor with N=None exposes "shaft") instead of an ordinary free
    parameter -- both are validated and closed the same way, just through
    free_parameters()/closes_parameters() for the former and
    free_mechanical_ports()/closes_mechanical_nodes() for the latter.

    Raises ValueError at construction if the target component doesn't
    currently declare free_param as free (in either sense) — this is almost
    always a configuration mistake (forgetting to leave that constructor arg
    as None), and failing fast here is far clearer than a mismatched Newton
    system discovered later as a solver error.
    """

    def __init__(
        self,
        name: str,
        component: BaseComponent,
        free_param: str,
        target_metric: str,
        value: TargetValue,
    ):
        mech_ports = component.mechanical_ports()
        self._mechanical = (
            free_param in mech_ports and mech_ports[free_param] in component.free_mechanical_ports()
        )
        if not self._mechanical and free_param not in component.free_parameters():
            raise ValueError(
                f"Setpoint {name!r} targets {component.name!r}.{free_param}, but that "
                f"component doesn't currently declare {free_param!r} as a free parameter "
                f"or a free mechanical port — pass None for it at construction so it "
                f"becomes a solvable unknown."
            )

        self.name = name
        self.component = component
        self.free_param = free_param
        self.target_metric = target_metric
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
        metrics = self.component.report_metrics(state)
        if metrics is None or self.target_metric not in metrics:
            raise ValueError(
                f"Setpoint {self.name!r} targets metric {self.target_metric!r} on "
                f"{self.component.name!r}, but report_metrics() doesn't expose it "
                f"(got: {sorted(metrics) if metrics else []})"
            )
        return [metrics[self.target_metric] - self._target(state)]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        measured = self.component.report_metrics(state)[self.target_metric]
        target = self._target(state)
        return {
            "target [-]": target,
            "measured [-]": measured,
            "error [-]": measured - target,
        }
