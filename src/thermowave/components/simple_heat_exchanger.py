from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_MIN_MDOT = 1.0e-9  # kg/s, floor for the Q/mdot energy-residual division below --
# a Newton iterate that clamps a free mdot near zero (Solver.MDOT_MIN) would
# otherwise blow this up.


class SimpleHeatExchanger(BaseComponent):
    """Single-stream heat addition/removal, 0D model: one fluid network, one
    duty Q.

    Two ports (in/out) -- this is not a two-stream exchanger (see
    HeatExchanger/MultiPassHeatExchanger for that); it's the "apply Q to
    this stream" building block, e.g. a heater/cooler/duty specified
    directly rather than derived from a second stream's own state.

    Q [W]: positive heats the fluid, negative cools it -- give it directly,
    or leave it None to make it a free Newton unknown (closed by a
    Setpoint/Controller/PIDController targeting some downstream quantity,
    the same free-parameter pattern as Combustor's mdot_fuel). Not clamped
    either way; a Q that drives the outlet temperature outside the fluid
    model's valid range fails the same way any other component's residuals
    would.

    Pressure drop is a simple fixed ratio (P_out = PR * P_in, same style as
    SimpleCompressor/SimpleTurbine), not a K-factor loss model.
    """

    def __init__(self, name: str, Q: float | None, PR: float = 1.0):
        if PR <= 0:
            raise ValueError(f"SimpleHeatExchanger {name!r}: PR must be > 0, got {PR}")
        self.name = name
        self.Q = Q
        self.PR = PR
        self._in_node = f"{name}.in"
        self._out_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"in": self._in_node, "out": self._out_node}

    def report_category(self) -> str:
        return "single_stream_hx"

    def free_parameters(self) -> dict[str, float]:
        if self.Q is not None:
            return {}
        return {"Q": 0.0}

    def _duty(self, state: "NetworkState") -> float:
        if self.Q is not None:
            return self.Q
        return state.param(f"{self.name}.Q")

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        Q_guess = self.Q if self.Q is not None else 0.0
        return self.PR * P_in, h_in + Q_guess / max(mdot, _MIN_MDOT)

    def residuals(self, state: "NetworkState") -> list[float]:
        Q = self._duty(state)
        P_in, h_in = state.node(self._in_node)
        P_out, h_out = state.node(self._out_node)
        mdot = state.mdot(self._in_node)

        momentum_residual = P_out - self.PR * P_in
        energy_residual = h_out - (h_in + Q / max(mdot, _MIN_MDOT))
        mass_residual = state.mdot(self._out_node) - mdot
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._in_node)
        P_out, h_out = state.node(self._out_node)
        return {
            "power [W]": self._duty(state),
            "PR [-]": self.PR,
            "T_in [K]": state.fluid_at(self._in_node).temperature_ph(P_in, h_in),
            "T_out [K]": state.fluid_at(self._out_node).temperature_ph(P_out, h_out),
        }
