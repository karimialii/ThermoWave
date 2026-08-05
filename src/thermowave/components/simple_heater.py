from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class SimpleHeater(BaseComponent):
    """Single-stream heater targeting a fixed outlet temperature, without
    modeling the second stream that supplies the heat -- the single-phase
    -safe counterpart to SimpleEvaporator's outlet-spec mode.

    SimpleEvaporator/SimpleCondenser target a state *relative to the
    saturation curve* (superheat/subcool/quality), which is exactly wrong
    for a fluid that stays supercritical across the whole heater (e.g. sCO2
    cycle recuperator->heater crossings well above CO2's ~73.8 bar critical
    pressure): there is no saturation temperature to be relative TO up
    there, and CoolProp raises rather than returning one. SimpleHeater
    sidesteps this by targeting an outlet temperature directly --
    h_out = fluid.enthalpy_pt(P_out, T_out) needs nothing from the
    saturation dome, so it works identically whether the fluid ends up
    sub- or supercritical. (For a fluid that DOES stay subcritical, prefer
    SimpleEvaporator/SimpleCondenser instead when a superheat/subcool/
    quality target is the natural way to say what's wanted -- SimpleHeater
    is the fallback for when that framing doesn't apply, not a replacement
    for it.)

    Duty is NOT computed from an effectiveness/UA calculation (no second
    stream is modeled at all) -- it comes from the outlet state:
    Q = mdot * (h_out - h_in) is a *reported* result, not a residual, same
    convention as SimpleEvaporator/SimpleCondenser.

    P_out = PR * P_in (PR default 1.0; give PR < 1 for a heater that also
    carries a pressure drop, e.g. an external heater/receiver whose own
    pressure drop is known but not separately modeled). 3 residuals:
    momentum, energy, mass.
    """

    def __init__(self, name: str, T_out: float, PR: float = 1.0):
        if not (0.0 < PR <= 1.0):
            raise ValueError(f"SimpleHeater {name!r}: PR must be in (0, 1], got {PR}")
        if T_out <= 0.0:
            raise ValueError(f"SimpleHeater {name!r}: T_out must be > 0 (Kelvin), got {T_out}")
        self.name = name
        self.PR = PR
        self.T_out = T_out
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        return self.PR * P_in, h_in

    def residuals(self, state: "NetworkState") -> list[float]:
        fluid = state.fluid_at(self._inlet_node)
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        mdot_in = state.mdot(self._inlet_node)

        h_out_target = fluid.enthalpy_pt(P_out, self.T_out)

        momentum_residual = P_out - self.PR * P_in
        energy_residual = h_out - h_out_target
        mass_residual = state.mdot(self._outlet_node) - mdot_in
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        mdot = state.mdot(self._inlet_node)
        fluid = state.fluid_at(self._outlet_node)
        return {
            "power [W]": mdot * (h_out - h_in),
            "PR [-]": P_out / P_in,
            "T_out [K]": fluid.temperature_ph(P_out, h_out),
        }
