from __future__ import annotations

import math
from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_RAD_PER_MIN_TO_RPM = 60.0 / (2.0 * math.pi)


class Shaft(BaseComponent):
    """A mechanical connection: two or more components sharing one physical
    shaft (rigid, no ports of its own — a mechanical, not flow, connection).

    Two modes, chosen by dynamic:

    dynamic=False (default): the classic steady-state coupling. Every
    connected component's own free shaft-speed parameter (usually "N",
    declared free the same way as for Setpoint/Controller — see
    Compressor/Turbine's N=None) is tied to the first component's, scaled by
    a gear ratio: components[i].<free_param> == gear_ratios[i-1] *
    components[0].<free_param>, one residual per follower (gear_ratios
    defaults to all 1.0). components[0]'s own speed is still just a free
    Newton unknown with no residual of its own here — something else (a
    Setpoint or Controller) must still pin it down, exactly as before.

    dynamic=True: this Shaft owns its own rotor speed as differential state
    (via differential_parameters()/state_derivative() — see
    BaseComponent's docstrings), integrated from net torque and inertia:
    d(N)/dt = (net_power / omega) / inertia. Every connected component
    (including components[0] — there's no asymmetric "reference" component
    in this mode) is tied to that shared shaft speed: components[i].
    <free_param> == gear_ratios[i] * shaft_N, one residual per component
    (gear_ratios has one entry per component here, not per follower).
    Network.solve() closes shaft_N via "derivative == 0" (the torque-balance
    speed); Network.solve_transient() integrates it forward instead — no
    Setpoint/Controller needed for the shaft's own speed either way. N0
    seeds the initial guess (steady state) or the t=0 value (transient,
    unless an explicit `initial` SolveResult is given to solve_transient()).
    inertia must be > 0 in this mode.

    inertia [kg*m^2]: the shaft's (and everything rigidly attached to it)
    rotational inertia. Only used (and required > 0) when dynamic=True; in
    the default steady-state-only mode it's stored but has no effect,
    purely so switching a network over to dynamic=True later doesn't
    require re-specifying it.

    Torque-only members: a listed component that does NOT declare free_param
    as free (e.g. a ShaftLoad — a generator's commanded power draw) is a
    pure torque contribution: it has no speed unknown of its own, so it gets
    no speed-tie residual and no gear_ratios entry, but its
    report_metrics()["power [W]"] still enters the power balance with its
    sign. This is how an electrical load enters the physics: in dynamic
    mode the steady closing equation "derivative == 0" then reads
    turbine power == compressor power + load power (all scaled by
    efficiency/signs), so the equilibrium speed genuinely depends on the
    demanded power — and a transient integrates the imbalance through the
    rotor inertia. gear_ratios/speed-tie bookkeeping counts only the
    speed-tied components; at least one component must be speed-tied
    (dynamic mode), or the first component must be plus at least one more
    (static mode, where the first speed-tied component is the reference).

    efficiency: mechanical transmission efficiency (bearing/windage losses
    folded into one number, 1.0 = lossless). Not enforced as a residual —
    each connected machine's own residuals() already closes its own energy
    balance on the fluid side, so there's no separate path for lost power
    to go if efficiency were mixed into the coupling equation itself.
    Instead it scales the net shaft power reported by this component (and,
    in dynamic mode, the torque driving the speed integration): pass signs
    (default all +1.0) to mark which components deliver power to the shaft
    (e.g. a turbine: +1.0) vs draw from it (e.g. a compressor: -1.0), and
    report_metrics()["power [W]"] becomes
    efficiency * sum(sign_i * component_i power).
    """

    def __init__(
        self,
        name: str,
        components: list[BaseComponent],
        gear_ratios: list[float] | None = None,
        signs: list[float] | None = None,
        efficiency: float = 1.0,
        inertia: float = 0.0,
        free_param: str = "N",
        dynamic: bool = False,
        N0: float = 50000.0,
    ):
        if len(components) < 2:
            raise ValueError(f"Shaft {name!r} needs at least 2 components, got {len(components)}")
        # Components declaring free_param carry a speed unknown to tie to the
        # shaft; the rest (e.g. a ShaftLoad) are torque-only members that
        # enter the power balance but have no speed of their own.
        speed_tied = [c for c in components if free_param in c.free_parameters()]
        min_tied = 1 if dynamic else 2
        if len(speed_tied) < min_tied:
            raise ValueError(
                f"Shaft {name!r} needs at least {min_tied} component(s) declaring "
                f"{free_param!r} as free (got {len(speed_tied)}) — a shaft made only of "
                f"torque-only members (e.g. ShaftLoad) has no speed unknown to couple; "
                f"pass None for {free_param!r} on the turbomachines so it becomes a "
                f"solvable unknown."
            )
        n_ratios = len(speed_tied) if dynamic else len(speed_tied) - 1
        if gear_ratios is None:
            gear_ratios = [1.0] * n_ratios
        if len(gear_ratios) != n_ratios:
            raise ValueError(
                f"Shaft {name!r}: gear_ratios must have one entry per "
                f"{'speed-tied component' if dynamic else 'speed-tied follower'} "
                f"({n_ratios}), got {len(gear_ratios)}"
            )
        if signs is None:
            signs = [1.0] * len(components)
        if len(signs) != len(components):
            raise ValueError(
                f"Shaft {name!r}: signs must have one entry per component "
                f"({len(components)}), got {len(signs)}"
            )
        if dynamic and inertia <= 0.0:
            raise ValueError(
                f"Shaft {name!r} has dynamic=True, which integrates its own speed from "
                f"net torque and inertia; inertia={inertia} must be > 0."
            )

        self.name = name
        self.components = components
        self._speed_tied = speed_tied
        self.gear_ratios = gear_ratios
        self.signs = signs
        self.efficiency = efficiency
        self.inertia = inertia
        self.free_param = free_param
        self.dynamic = dynamic
        self.N0 = N0

    def ports(self) -> dict[str, str]:
        return {}

    def report_category(self) -> str:
        return "shaft"

    def _speed(self, state: "NetworkState", component: BaseComponent) -> float:
        return state.param(f"{component.name}.{self.free_param}")

    def _reference_speed(self, state: "NetworkState") -> float:
        if self.dynamic:
            return state.param(f"{self.name}.N")
        return self._speed(state, self._speed_tied[0])

    def _net_power(self, state: "NetworkState") -> float:
        net_power = 0.0
        for component, sign in zip(self.components, self.signs):
            metrics = component.report_metrics(state)
            if metrics is not None and "power [W]" in metrics:
                net_power += sign * metrics["power [W]"]
        return net_power * self.efficiency

    def differential_parameters(self) -> dict[str, float]:
        if not self.dynamic:
            return {}
        return {"N": self.N0}

    def state_derivative(self, state: "NetworkState") -> dict[str, float]:
        if not self.dynamic:
            return {}
        N = state.param(f"{self.name}.N")
        omega = N / _RAD_PER_MIN_TO_RPM
        net_power = self._net_power(state)
        torque = net_power / omega if omega > 0.0 else 0.0
        return {"N": (torque / self.inertia) * _RAD_PER_MIN_TO_RPM}

    def residuals(self, state: "NetworkState") -> list[float]:
        if self.dynamic:
            N_ref = self._reference_speed(state)
            return [
                self._speed(state, component) - ratio * N_ref
                for component, ratio in zip(self._speed_tied, self.gear_ratios)
            ]
        N_ref = self._speed(state, self._speed_tied[0])
        return [
            self._speed(state, component) - ratio * N_ref
            for component, ratio in zip(self._speed_tied[1:], self.gear_ratios)
        ]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        return {
            "N [rev/min]": self._reference_speed(state),
            "power [W]": self._net_power(state),
            "eta [-]": self.efficiency,
            "inertia [kg*m^2]": self.inertia,
        }
