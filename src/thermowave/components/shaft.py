from __future__ import annotations

import math
import warnings
from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_RAD_PER_MIN_TO_RPM = 60.0 / (2.0 * math.pi)


class Shaft(BaseComponent):
    """A mechanical connection: two or more components sharing one physical
    shaft, wired explicitly via Network.connect(..., kind="mechanical")/
    kind="signal" — this component itself contributes no flow ports.

    `members` lists every component riding on this shaft, in the order their
    ports are exposed: this Shaft gets one mechanical port ("m0", "m1", ...)
    per member that has its own "shaft" mechanical port (a real speed
    unknown — e.g. Turbine/Compressor, whether N is given or left free), and
    one signal port ("p0", "p1", ...) per member in `members` order (every
    member, including a torque-only one like ShaftLoad, since power still
    needs to reach the balance). The caller wires both explicitly after
    add_component():

        shaft = Shaft("sh1", members=[turb, comp, load], signs=[1.0, -1.0, -1.0])
        network.add_component(shaft)
        network.connect(shaft, "m0", turb, "shaft", kind="mechanical")
        network.connect(shaft, "m1", comp, "shaft", kind="mechanical")
        network.connect(shaft, "p0", turb, "power", kind="signal")
        network.connect(shaft, "p1", comp, "power", kind="signal")
        network.connect(shaft, "p2", load, "power", kind="signal")

    ShaftLoad has no "shaft" mechanical port (no speed unknown of its own —
    a pure torque contribution), so it gets no mechanical port here, only a
    signal one.

    Two modes, chosen by dynamic:

    dynamic=False (default): the classic steady-state coupling. Every
    speed-tied member's shaft speed is tied to the first speed-tied member's
    ("m0"), scaled by a gear ratio: N("m{i}") == gear_ratios[i-1] * N("m0"),
    one residual per follower (gear_ratios defaults to all 1.0). "m0"'s own
    speed is still just whatever Turbine/Compressor/connect() made it (fixed
    or a free Newton unknown) — something else (a Setpoint/Controller, or a
    fixed N on that member) must still pin it down if it's free, exactly as
    before.

    dynamic=True: this Shaft owns its own rotor speed as differential state
    (via differential_parameters()/state_derivative()), integrated from net
    torque and inertia: d(N)/dt = (net_power / omega) / inertia. Every
    speed-tied member (there's no asymmetric "reference" member in this
    mode) is tied to that shared shaft speed: N("m{i}") == gear_ratios[i] *
    shaft_N, one residual per member (gear_ratios has one entry per
    speed-tied member here, not per follower). Network.solve() closes
    shaft_N via "derivative == 0" (the torque-balance speed);
    Network.solve_transient() integrates it forward instead — no
    Setpoint/Controller needed for the shaft's own speed either way. N0
    seeds the initial guess (steady state) or the t=0 value (transient,
    unless an explicit `initial` SolveResult is given to solve_transient()).
    inertia must be > 0 in this mode.

    inertia [kg*m^2]: the shaft's (and everything rigidly attached to it)
    rotational inertia. Only used (and required > 0) when dynamic=True; in
    the default steady-state-only mode it's stored but has no effect,
    purely so switching a network over to dynamic=True later doesn't
    require re-specifying it.

    Torque-only members (e.g. ShaftLoad — a generator's commanded power
    draw): no speed unknown of its own, so it gets no mechanical port, no
    speed-tie residual, and no gear_ratios entry, but its published "power"
    signal still enters the power balance with its sign. This is how an
    electrical load enters the physics: in dynamic mode the steady closing
    equation "derivative == 0" then reads turbine power == compressor power
    + load power (all scaled by efficiency/signs), so the equilibrium speed
    genuinely depends on the demanded power — and a transient integrates the
    imbalance through the rotor inertia. gear_ratios/speed-tie bookkeeping
    counts only the speed-tied members; at least one member must be
    speed-tied (dynamic mode), or the first speed-tied member plus at least
    one more (static mode, where the first speed-tied member is the
    reference).

    efficiency: mechanical transmission efficiency (bearing/windage losses
    folded into one number, 1.0 = lossless). Not enforced as a residual —
    each connected machine's own residuals() already closes its own energy
    balance on the fluid side, so there's no separate path for lost power
    to go if efficiency were mixed into the coupling equation itself.
    Instead it scales the net shaft power reported by this component (and,
    in dynamic mode, the torque driving the speed integration): pass signs
    (default all +1.0) to mark which members deliver power to the shaft
    (e.g. a turbine: +1.0) vs draw from it (e.g. a compressor: -1.0), and
    report_metrics()["power [W]"] becomes
    efficiency * sum(sign_i * member_i power).
    """

    def __init__(
        self,
        name: str,
        members: list[BaseComponent],
        gear_ratios: list[float] | None = None,
        signs: list[float] | None = None,
        efficiency: float = 1.0,
        inertia: float = 0.0,
        dynamic: bool = False,
        N0: float = 50000.0,
    ):
        if len(members) < 2:
            raise ValueError(f"Shaft {name!r} needs at least 2 members, got {len(members)}")
        # Members with a "shaft" mechanical port carry a real speed unknown
        # to tie together; the rest (e.g. ShaftLoad) are torque-only and
        # enter only the power balance, never the speed ties.
        speed_tied = [c for c in members if "shaft" in c.mechanical_ports()]
        min_tied = 1 if dynamic else 2
        if len(speed_tied) < min_tied:
            raise ValueError(
                f"Shaft {name!r} needs at least {min_tied} member(s) with their own "
                f"'shaft' mechanical port (got {len(speed_tied)}) — a shaft made only of "
                f"torque-only members (e.g. ShaftLoad) has no speed unknown to couple."
            )
        n_ratios = len(speed_tied) if dynamic else len(speed_tied) - 1
        if gear_ratios is None:
            gear_ratios = [1.0] * n_ratios
        if len(gear_ratios) != n_ratios:
            raise ValueError(
                f"Shaft {name!r}: gear_ratios must have one entry per "
                f"{'speed-tied member' if dynamic else 'speed-tied follower'} "
                f"({n_ratios}), got {len(gear_ratios)}"
            )
        if signs is None:
            signs = [1.0] * len(members)
        if len(signs) != len(members):
            raise ValueError(
                f"Shaft {name!r}: signs must have one entry per member "
                f"({len(members)}), got {len(signs)}"
            )
        if dynamic and inertia <= 0.0:
            raise ValueError(
                f"Shaft {name!r} has dynamic=True, which integrates its own speed from "
                f"net torque and inertia; inertia={inertia} must be > 0."
            )
        if dynamic and efficiency != 1.0:
            warnings.warn(
                f"Shaft {name!r} sets efficiency={efficiency} with dynamic=True, where it "
                f"has almost no effect: efficiency scales the *net* shaft power, and the "
                f"steady equilibrium a dynamic shaft settles at is exactly net == 0, so "
                f"the factor cancels and the mechanical loss silently disappears. Put it "
                f"on the load instead — ShaftLoad(efficiency=eta_gen * eta_mech) with "
                f"Shaft(efficiency=1.0) — so the shaft-side draw becomes "
                f"P_elec / (eta_gen * eta_mech) and the loss is actually carried.",
                stacklevel=2,
            )

        self.name = name
        self.members = members
        self._speed_tied = speed_tied
        self.gear_ratios = gear_ratios
        self.signs = signs
        self.efficiency = efficiency
        self.inertia = inertia
        self.dynamic = dynamic
        self.N0 = N0
        self._mech_ports = {f"m{i}": f"{name}.m{i}" for i in range(len(speed_tied))}
        self._signal_port_names = {f"p{i}": f"{name}.p{i}" for i in range(len(members))}

    def ports(self) -> dict[str, str]:
        return {}

    def mechanical_ports(self) -> dict[str, str]:
        return dict(self._mech_ports)

    def signal_ports(self) -> dict[str, str]:
        return dict(self._signal_port_names)

    def closes_mechanical_nodes(self) -> list[str]:
        """The canonical mechanical nodes this shaft's speed-tie residuals
        pin down.

        dynamic=True owns its own speed as differential state and ties every
        speed-tied member to it, so all of them are closed here. dynamic=False
        ties followers to "m0"'s speed instead, which leaves that reference
        speed itself for something external (a Setpoint/Controller, or a
        fixed N on that member) to close — so "m0" is deliberately absent.
        """
        network = getattr(self, "_network", None)
        if network is None:
            return []
        names = (
            list(self._mech_ports) if self.dynamic else list(self._mech_ports)[1:]
        )
        return [network._canonical(self._mech_ports[n]) for n in names]

    def report_category(self) -> str:
        return "shaft"

    def _reference_speed(self, state: "NetworkState") -> float:
        if self.dynamic:
            return state.param(f"{self.name}.N")
        return state.N(self._mech_ports["m0"])

    def _net_power(self, state: "NetworkState") -> float:
        net_power = 0.0
        for i, sign in enumerate(self.signs):
            port_id = self._signal_port_names[f"p{i}"]
            net_power += sign * state.signal(port_id)
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
        mech_ids = list(self._mech_ports.values())
        if self.dynamic:
            N_ref = self._reference_speed(state)
            return [
                state.N(port_id) - ratio * N_ref
                for port_id, ratio in zip(mech_ids, self.gear_ratios)
            ]
        N_ref = state.N(mech_ids[0])
        return [
            state.N(port_id) - ratio * N_ref
            for port_id, ratio in zip(mech_ids[1:], self.gear_ratios)
        ]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        return {
            "N [rev/min]": self._reference_speed(state),
            "power [W]": self._net_power(state),
            "eta [-]": self.efficiency,
            "inertia [kg*m^2]": self.inertia,
        }
