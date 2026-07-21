from __future__ import annotations

from typing import TYPE_CHECKING, Union

from thermowave.components.base_component import BaseComponent
from thermowave.core.constants import STEFAN_BOLTZMANN

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

# One endpoint of a heat path: a ThermalMass (its own differential
# temperature), a fixed float (e.g. ambient), or (component, port_name)
# reading that component's live fluid temperature at that port's node.
TemperatureSource = Union["ThermalMass", float, tuple[BaseComponent, str]]


def _temperature_of(source: TemperatureSource, state: "NetworkState") -> float:
    if isinstance(source, ThermalMass):
        return state.param(f"{source.name}.T")
    if isinstance(source, (int, float)):
        return float(source)
    component, port = source
    node = component.ports()[port]
    return state.fluid_at(node).temperature_ph(*state.node(node))


def heat_loss_watts(heat_path: BaseComponent | None, state: "NetworkState") -> float:
    """Q [W] a flow component (Turbine/Compressor/Combustor/...) is currently
    losing through its optional heat_path (Convection/Conduction/Radiation),
    or 0.0 if it has none. Positive = heat leaving that component's fluid,
    same sign convention Pipe's own heat_loss already uses."""
    if heat_path is None:
        return 0.0
    return heat_path.Q(state)


class ThermalMass(BaseComponent):
    """A solid's own temperature as time-integrated state (a casing, a
    shaft, ...) — the node any Convection/Conduction/Radiation heat path
    attaches to on the solid side. Has no flow ports (ports() -> {}, like
    Generator): it never joins the fluid network graph, only the Newton
    system's differential-state bookkeeping.

    thermal_capacitance is the lumped m*cp [J/K] of the solid, given
    directly rather than derived from a separate mass/material-cp pair —
    consistent with how MultiPassHeatExchanger's UA is already an opaque
    input here rather than computed from geometry.

    heat_sources is a plain mutable list of (heat_path, sign) pairs, left
    empty at construction and appended to (or assigned wholesale) once the
    Convection/Conduction/Radiation components referencing this mass exist
    — they must be built after both of their endpoints, so this mass can't
    know about them yet at its own construction time. sign is +1.0 if this
    mass is heat_path's `b` endpoint (gains heat_path.Q()) or -1.0 if it's
    the `a` endpoint (loses it) — same convention Shaft's components/signs
    pair already uses for summing signed contributions from other
    components' report_metrics().

    dT/dt = (sum of sign*heat_path.Q(state) over heat_sources) /
    thermal_capacitance. Network.solve() closes T to whatever value makes
    that sum zero (steady state); Network.solve_transient() integrates T
    forward from T0 instead.
    """

    def __init__(self, name: str, thermal_capacitance: float, T0: float):
        if thermal_capacitance <= 0.0:
            raise ValueError(
                f"ThermalMass {name!r}: thermal_capacitance must be > 0, "
                f"got {thermal_capacitance}"
            )
        self.name = name
        self.thermal_capacitance = thermal_capacitance
        self.T0 = T0
        self.heat_sources: list[tuple[BaseComponent, float]] = []

    def ports(self) -> dict[str, str]:
        return {}

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def report_category(self) -> str:
        return "thermal_mass"

    def differential_parameters(self) -> dict[str, float]:
        return {"T": self.T0}

    def _net_heat(self, state: "NetworkState") -> float:
        return sum(sign * heat_path.Q(state) for heat_path, sign in self.heat_sources)

    def state_derivative(self, state: "NetworkState") -> dict[str, float]:
        return {"T": self._net_heat(state) / self.thermal_capacitance}

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        return {
            "T [K]": state.param(f"{self.name}.T"),
            "Q_net [W]": self._net_heat(state),
        }


def _heat_path_metrics(path: "Convection | Conduction | Radiation", state: "NetworkState") -> dict[str, float]:
    return {
        "Q [W]": path.Q(state),
        "T_a [K]": _temperature_of(path.a, state),
        "T_b [K]": _temperature_of(path.b, state),
    }


class Convection(BaseComponent):
    """Q = h*A*(T_a - T_b): the linear convective heat-transfer law, for
    either free or forced convection — they differ only in how h is
    physically obtained (natural buoyancy-driven flow vs. an imposed flow
    over a surface), not in this formula. h is given directly here (no
    Nusselt-Rayleigh/Reynolds correlation yet); a future helper can derive
    h from geometry and hand it to this same class unchanged.

    No flow ports (ports() -> {}, like Generator) and no algebraic
    residuals — it's a passive reader other components' state_derivative()/
    residuals() pull Q(state) from, not a network graph node. Q is positive
    when a is hotter than b (heat flowing from a to b), the same convention
    Pipe's own heat_loss already uses ("positive = heat lost" by whichever
    side is endpoint a).
    """

    def __init__(self, name: str, a: TemperatureSource, b: TemperatureSource, h: float, A: float):
        self.name = name
        self.a = a
        self.b = b
        self.h = h
        self.A = A

    def ports(self) -> dict[str, str]:
        return {}

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def report_category(self) -> str:
        return "heat_transfer"

    def Q(self, state: "NetworkState") -> float:
        return self.h * self.A * (_temperature_of(self.a, state) - _temperature_of(self.b, state))

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        return _heat_path_metrics(self, state)


class Conduction(BaseComponent):
    """Q = (k*A/L)*(T_a - T_b): steady 1D conduction through a solid path
    of conductivity k, cross-sectional area A, and length L (e.g. a shaft
    conducting heat from a turbine casing to a compressor casing). Same
    zero-port, zero-residual, a->b sign convention as Convection."""

    def __init__(
        self, name: str, a: TemperatureSource, b: TemperatureSource, k: float, A: float, L: float
    ):
        self.name = name
        self.a = a
        self.b = b
        self.k = k
        self.A = A
        self.L = L
        self.UA = k * A / L

    def ports(self) -> dict[str, str]:
        return {}

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def report_category(self) -> str:
        return "heat_transfer"

    def Q(self, state: "NetworkState") -> float:
        return self.UA * (_temperature_of(self.a, state) - _temperature_of(self.b, state))

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        return _heat_path_metrics(self, state)


class Radiation(BaseComponent):
    """Q = emissivity*view_factor*sigma*A*(T_a^4 - T_b^4): surface-to-
    surface (or surface-to-ambient) radiative exchange. General-purpose,
    not combustor-specific — the intended primitive a future 1D combustion-
    chamber liner model builds on, discretizing a liner into many of these
    rather than needing its own radiation formula. Same zero-port,
    zero-residual, a->b sign convention as Convection/Conduction.
    """

    def __init__(
        self,
        name: str,
        a: TemperatureSource,
        b: TemperatureSource,
        emissivity: float,
        A: float,
        view_factor: float = 1.0,
    ):
        if not (0.0 < emissivity <= 1.0):
            raise ValueError(
                f"Radiation {name!r}: emissivity must be in (0, 1], got {emissivity}"
            )
        self.name = name
        self.a = a
        self.b = b
        self.emissivity = emissivity
        self.A = A
        self.view_factor = view_factor

    def ports(self) -> dict[str, str]:
        return {}

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def report_category(self) -> str:
        return "heat_transfer"

    def Q(self, state: "NetworkState") -> float:
        T_a = _temperature_of(self.a, state)
        T_b = _temperature_of(self.b, state)
        return (
            self.emissivity * self.view_factor * STEFAN_BOLTZMANN * self.A
            * (T_a**4 - T_b**4)
        )

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        return _heat_path_metrics(self, state)
