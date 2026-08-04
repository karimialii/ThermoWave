from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.components.sensor import SENSOR_QUANTITIES as _ALLOWED_QUANTITIES

if TYPE_CHECKING:
    from collections.abc import Callable

    from thermowave.core.network import NetworkState


class PIDController(BaseComponent):
    """Time-domain PID controller, for use inside Network.solve_transient().

    Controller pins its target exactly on every steady-state solve — an
    ideal, infinite-gain controller with no dynamics of its own. PIDController
    is the finite-response counterpart: it holds free_param at a fixed
    self.output during each algebraic solve (one more Newton residual pinning
    it, same shape as Controller's), but that output only *changes* once per
    transient step, via step() computing Kp*error + Ki*integral +
    Kd*derivative from the current Sensor reading — not solved for exactly
    on every call. That's what lets it show overshoot, offset, and finite
    settling time the way a real controller does.

    step() is not part of the Newton system — it's an explicit-time update,
    called once per step by Network.solve_transient() (the same role Shaft's
    rotor-speed integration plays there), using the state from the step that
    was just solved to compute the output for the *next* one.

    error = setpoint - measured; sign/magnitude of Kp determines whether
    increasing free_param increases or decreases the measured quantity —
    get this backwards and the loop diverges instead of converging (same
    "check the pairing makes physical sense" caveat as tuning a real PID).

    output0 both seeds self.output for the very first solve (must already put
    free_param somewhere reasonable) and acts as a fixed bias added to every
    later output: output = bias + Kp*error + Ki*integral + Kd*derivative.
    Without that bias, a free_param whose sensible operating range doesn't
    straddle zero (a shaft speed in rev/min, say) could only be reached by
    winding the integral term up to that whole magnitude, which fights
    anti-windup clamping and makes startup slow/oscillatory; the bias lets
    Kp/Ki/Kd instead represent a small correction around a known-reasonable
    starting point. output_min/output_max clamp self.output if given; the
    integral term uses simple clamped (anti-windup) integration — it only
    accumulates further in the direction that's already saturating once
    output_min/output_max is hit, so it doesn't keep growing unboundedly
    while saturated.

    feedforward: an optional callable(state) -> float replacing that constant
    bias with one computed fresh each step, so the operating point the PID
    trims around can *move* with the plant instead of being pinned to
    whatever it was at t=0. This matters wherever a machine's required
    actuator setting is mostly a known function of its load: on a gas turbine
    the fuel flow needed is set by the energy balance, so biasing from the
    commanded load (feedforward=lambda state: fuel_for(load.power)) leaves
    the PI trimming a small correction, rather than having to integrate
    across the whole actuator range to find the operating point at all. The
    difference isn't just settling time — a plant can have regions with no
    steady state at all (too little fuel to carry the load), and an integral
    term wandering to find its bias can walk straight into one. self.bias is
    used when feedforward is None.
    """

    def __init__(
        self,
        name: str,
        sensor: BaseComponent,
        quantity: str,
        component: BaseComponent,
        free_param: str,
        setpoint: float,
        Kp: float,
        Ki: float = 0.0,
        Kd: float = 0.0,
        output0: float = 0.0,
        output_min: float | None = None,
        output_max: float | None = None,
        feedforward: "Callable[[NetworkState], float] | None" = None,
    ):
        mech_ports = component.mechanical_ports()
        self._mechanical = (
            free_param in mech_ports and mech_ports[free_param] in component.free_mechanical_ports()
        )
        if not self._mechanical and free_param not in component.free_parameters():
            raise ValueError(
                f"PIDController {name!r} actuates {component.name!r}.{free_param}, but "
                f"that component doesn't currently declare {free_param!r} as a free "
                f"parameter or a free mechanical port — pass None for it at "
                f"construction so it becomes a solvable unknown."
            )

        self.name = name
        self.sensor = sensor
        self.quantity = quantity
        self.component = component
        self.free_param = free_param
        self.setpoint = setpoint
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.output_min = output_min
        self.output_max = output_max
        self.feedforward = feedforward
        self.output = output0
        # Public so a Schedule can drive it along a time profile
        # (Schedule(target=pid, attr="bias", ...)) the same way it drives a
        # setpoint — the simplest form of feedforward when the operating
        # point is known as a function of time rather than of state.
        self.bias = output0
        self._output0 = output0
        self._integral = 0.0
        self._prev_error: float | None = None

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

    def _clamp(self, value: float) -> float:
        if self.output_min is not None:
            value = max(self.output_min, value)
        if self.output_max is not None:
            value = min(self.output_max, value)
        return value

    def _measured(self, state: "NetworkState") -> float:
        metrics = self.sensor.report_metrics(state)
        if self.quantity not in metrics:
            raise ValueError(
                f"PIDController {self.name!r} reads quantity {self.quantity!r} from "
                f"sensor {self.sensor.name!r}, but it doesn't expose that reading "
                f"(got: {sorted(metrics)}; try one of {_ALLOWED_QUANTITIES})"
            )
        return metrics[self.quantity]

    def reset(self, output0: float | None = None) -> None:
        """Clear the integral term, the remembered previous error, and the
        current output — everything step() accumulates.

        A PIDController carries that state on the instance, so reusing one
        across two solve_transient() runs would otherwise silently start the
        second run with the first one's wound-up integral. output0 re-seeds
        both output and bias if given; otherwise both return to the value
        passed at construction.
        """
        if output0 is not None:
            self._output0 = output0
        self.output = self._output0
        self.bias = self._output0
        self._integral = 0.0
        self._prev_error = None

    def step(self, state: "NetworkState", dt: float) -> float:
        """Advance the controller by one transient step: read the sensor off
        the just-solved state, compute the PID law, and update self.output
        (the value residuals() will pin free_param to on the *next* solve).
        Returns the new self.output. dt <= 0 skips the derivative term (only
        relevant for the very first call, where there's no previous error).
        """
        error = self.setpoint - self._measured(state)
        derivative = 0.0
        if self._prev_error is not None and dt > 0.0:
            derivative = (error - self._prev_error) / dt
        candidate_integral = self._integral + error * dt if dt > 0.0 else self._integral

        bias = self.feedforward(state) if self.feedforward is not None else self.bias
        raw_output = (
            bias + self.Kp * error + self.Ki * candidate_integral + self.Kd * derivative
        )
        clamped_output = self._clamp(raw_output)
        if clamped_output == raw_output:
            # Not saturated: keep the integral update. Saturated: discard it
            # (anti-windup) so the integral term doesn't keep growing while
            # the output can't actually move any further in that direction.
            self._integral = candidate_integral

        self.output = clamped_output
        self._prev_error = error
        return self.output

    def residuals(self, state: "NetworkState") -> list[float]:
        if self._mechanical:
            port_id = self.component.mechanical_ports()[self.free_param]
            actual = state.N(port_id)
        else:
            actual = state.param(f"{self.component.name}.{self.free_param}")
        return [actual - self.output]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        measured = self._measured(state)
        return {
            "target [-]": self.setpoint,
            "measured [-]": measured,
            "error [-]": self.setpoint - measured,
            "output [-]": self.output,
        }
