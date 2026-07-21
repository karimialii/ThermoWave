from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class Setpoint(BaseComponent):
    """Drives one of a component's report_metrics() outputs to a target value
    by leaning on a free parameter that component already declared.

    Has no ports of its own — it contributes exactly one residual:
    component.report_metrics(state)[target_metric] - value == 0. It doesn't
    "free" anything itself; the target component must already expose the
    named parameter as free (e.g. Compressor(..., N=None)), which is where
    the matching extra Newton unknown comes from. Setpoint just supplies the
    other half: the equation that pins that unknown down.

    This generalizes "give a target instead of a direct input" to any
    component/metric pair without teaching every component a bespoke
    N-or-PR-or-power constructor: e.g. tie a Compressor's free N to a target
    power, tie it to a target PR, or (once other components expose their own
    free_parameters()) tie a Valve's opening to a target downstream pressure.

    Raises ValueError at construction if the target component doesn't
    currently declare free_param as free — this is almost always a
    configuration mistake (forgetting to leave that constructor arg as
    None), and failing fast here is far clearer than a mismatched Newton
    system discovered later as a solver error.
    """

    def __init__(
        self,
        name: str,
        component: BaseComponent,
        free_param: str,
        target_metric: str,
        value: float,
    ):
        if free_param not in component.free_parameters():
            raise ValueError(
                f"Setpoint {name!r} targets {component.name!r}.{free_param}, but that "
                f"component doesn't currently declare {free_param!r} as free — pass "
                f"None for it at construction so it becomes a solvable unknown."
            )

        self.name = name
        self.component = component
        self.free_param = free_param
        self.target_metric = target_metric
        self.value = value

    def ports(self) -> dict[str, str]:
        return {}

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
        return [metrics[self.target_metric] - self.value]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        measured = self.component.report_metrics(state)[self.target_metric]
        return {
            "target [-]": self.value,
            "measured [-]": measured,
            "error [-]": measured - self.value,
        }
