from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class Schedule(BaseComponent):
    """Drives a plain attribute on some other object (typically a
    Controller/PIDController/Setpoint's setpoint/value) along a piecewise
    time profile, during Network.solve_transient().

    Contributes no ports and no residuals — it isn't part of the Newton
    system at all. Instead it exposes step(state, dt), the same hook
    PIDController uses: solve_transient() discovers and calls it once per
    timestep (see thermowave.core.transient.solve_transient()), after that
    step's algebraic solve, to update the target attribute for the *next*
    step. At construction it also writes the t=0 value immediately, so the
    very first solve (before any step() call) already sees the right value.

    This replaces the common tutorial pattern of a hand-written Python loop
    around solve_transient() that manually reassigns
    `some_controller.setpoint = profile(t)` every step — the profile lives
    in the network instead, the same way a Shaft(dynamic=True) owns its own
    time integration instead of the caller stepping it by hand.

    breakpoints is a list of (t, value) pairs, sorted by t (strictly
    increasing — no repeated times). Before the first breakpoint or after
    the last, the value holds flat at the nearest endpoint. Between
    breakpoints: linear interpolation (interpolation="linear", the default)
    or a step held at the earlier breakpoint's value until the next one is
    reached (interpolation="step") — e.g. for a setpoint that should jump
    rather than ramp.
    """

    def __init__(
        self,
        name: str,
        target: object,
        attr: str,
        breakpoints: list[tuple[float, float]],
        interpolation: str = "linear",
    ):
        if len(breakpoints) < 2:
            raise ValueError(f"Schedule {name!r} needs at least two (t, value) breakpoints")
        times = [t for t, _ in breakpoints]
        if times != sorted(times) or len(set(times)) != len(times):
            raise ValueError(
                f"Schedule {name!r} breakpoints must be sorted by strictly increasing time"
            )
        if interpolation not in ("linear", "step"):
            raise ValueError(
                f"Schedule {name!r} interpolation must be 'linear' or 'step', "
                f"got {interpolation!r}"
            )
        if not hasattr(target, attr):
            raise ValueError(
                f"Schedule {name!r} targets {target!r}.{attr}, but that attribute "
                f"doesn't exist yet — set it once at construction on the target first."
            )

        self.name = name
        self.target = target
        self.attr = attr
        self.breakpoints = list(breakpoints)
        self.interpolation = interpolation
        self._t = 0.0
        setattr(self.target, self.attr, self.value_at(self._t))

    def ports(self) -> dict[str, str]:
        return {}

    def report_category(self) -> str:
        return "controller"

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def value_at(self, t: float) -> float:
        """Read-only lookup of what this schedule's value is (or will be) at
        time t, independent of self._t / step() — doesn't touch the target
        attribute. Useful for plotting/reporting the commanded profile
        itself (e.g. alongside a PIDController's actual tracked response),
        since the profile is a pure function of t but the live target
        attribute only ever holds its current value, not its history."""
        breakpoints = self.breakpoints
        if t <= breakpoints[0][0]:
            return breakpoints[0][1]
        if t >= breakpoints[-1][0]:
            return breakpoints[-1][1]
        for (t0, v0), (t1, v1) in zip(breakpoints, breakpoints[1:]):
            if t0 <= t < t1:
                if self.interpolation == "step":
                    return v0
                return v0 + (t - t0) / (t1 - t0) * (v1 - v0)
        return breakpoints[-1][1]  # unreachable given the checks above

    def step(self, state: "NetworkState", dt: float) -> float:
        """Advance the schedule by one transient step and write the new
        value onto the target attribute — mirrors PIDController.step()'s
        role (called once per timestep by solve_transient(), using elapsed
        dt rather than an absolute clock passed in from outside)."""
        self._t += dt
        value = self.value_at(self._t)
        setattr(self.target, self.attr, value)
        return value

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        return {"target [-]": getattr(self.target, self.attr)}
