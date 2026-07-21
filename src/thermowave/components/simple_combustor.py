from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.components.heat_transfer import heat_loss_watts
from thermowave.core.constants import MDOT_FUEL_GUESS_FRACTION

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid


class SimpleCombustor(BaseComponent):
    """Combustor with a fixed lower heating value (LHV) heat-release model.

    No combustion chemistry — fuel is treated as pure heat release into the
    existing working fluid (same single-fluid-per-network simplification as
    the rest of this codebase; see Combustor for a chemistry-based
    alternative using Cantera equilibrium products). Fuel mass is added to
    the flow (outlet mdot = inlet mdot + mdot_fuel):
        Q = mdot_fuel * LHV * efficiency
        mdot_out * h_out = mdot_in * h_in [+ mdot_fuel * h_fuel_in] + Q
        P_out = PR * P_in           (PR <= 1, combustor pressure loss)

    Two ways to get mdot_fuel, chosen by use_fuel_port:

    - use_fuel_port=False (default): mdot_fuel [kg/s] is given directly, or
      left None to drive the combustor by some other known quantity instead
      (a target turbine inlet temperature, a target heat release, ...) — it
      then becomes an extra Newton unknown (via free_parameters(), seeded
      from a generic ~2% fuel-air-ratio guess) and needs a matching residual
      from somewhere else in the network, e.g. a Setpoint or Controller
      tying report_metrics()["power [W]"] or a downstream Sensor's
      temperature reading to a target — same pattern as Compressor/
      Turbine's free N. Fuel's own sensible enthalpy is ignored (assumed
      negligible next to LHV, the standard simplification).
    - use_fuel_port=True: a genuine third port ("fuel_in") is added —
      connect a real fuel-supply branch to it (e.g. Source -> Pipe ->
      combustor, representing an actual fuel line premixing with the air
      just before combustion). mdot_fuel is then read directly from that
      port's own solved mdot (state.mdot(fuel_in)) — it's just another flow
      branch in the network, not a component-owned free parameter — and its
      actual (P, h) is used for a more complete energy balance (fuel's own
      sensible enthalpy is included, not just its LHV). mdot_fuel (the
      constructor arg) is ignored in this mode. Note: since this network
      still tracks state through one shared BaseFluid, the fuel stream is
      given the same fluid properties as the air stream — there's no
      separate fuel-gas property model, same limitation as everywhere else
      in this codebase that mixes streams.

    heat_path: optional Convection/Conduction/Radiation
    (thermowave.components.heat_transfer) representing heat lost through
    this combustor's liner/casing to something else. None (the default)
    means fully adiabatic (apart from efficiency < 1, which already models
    a form of loss — see the docstring above), unchanged from before this
    existed. Unlike Turbine/Compressor's per-unit-mass energy residual,
    this component's energy_residual is already power-valued (Q above is a
    total wattage, not a specific enthalpy), so heat_path's Q(state) is
    subtracted directly here with no /mdot division. Since the path needs
    (self, "out") as one of its own endpoints, it can only be built after
    this SimpleCombustor already exists — pass it here if you have it, or
    just set combustor.heat_path = path afterwards.
    """

    def __init__(
        self,
        name: str,
        LHV: float,
        PR: float = 0.97,
        efficiency: float = 1.0,
        mdot_fuel: float | None = None,
        use_fuel_port: bool = False,
        heat_path: BaseComponent | None = None,
    ):
        if not (0.0 < PR <= 1.0):
            raise ValueError(f"SimpleCombustor {name!r}: PR must be in (0, 1], got {PR}")
        if not (0.0 < efficiency <= 1.0):
            raise ValueError(
                f"SimpleCombustor {name!r}: efficiency must be in (0, 1], got {efficiency}"
            )
        self.name = name
        self.LHV = LHV
        self.PR = PR
        self.efficiency = efficiency
        self.mdot_fuel = mdot_fuel
        self.use_fuel_port = use_fuel_port
        self.heat_path = heat_path
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"
        self._fuel_in_node = f"{name}.fuel_in"

    def ports(self) -> dict[str, str]:
        ports = {"in": self._inlet_node, "out": self._outlet_node}
        if self.use_fuel_port:
            ports["fuel_in"] = self._fuel_in_node
        return ports

    def report_category(self) -> str:
        return "combustor"

    def _fuel_flow(self, state: "NetworkState") -> float:
        if self.use_fuel_port:
            return state.mdot(self._fuel_in_node)
        if self.mdot_fuel is not None:
            return self.mdot_fuel
        return state.param(f"{self.name}.mdot_fuel")

    def free_parameters(self) -> dict[str, float]:
        if self.use_fuel_port or self.mdot_fuel is not None:
            return {}
        return {"mdot_fuel": MDOT_FUEL_GUESS_FRACTION}

    def guess_free_parameters(
        self, fluid: "BaseFluid", P_in: float, h_in: float, mdot: float
    ) -> dict[str, float]:
        if self.use_fuel_port or self.mdot_fuel is not None:
            return {}
        return {"mdot_fuel": MDOT_FUEL_GUESS_FRACTION * mdot}

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        fuel_guess = self.mdot_fuel if self.mdot_fuel is not None else MDOT_FUEL_GUESS_FRACTION * mdot
        Q_guess = fuel_guess * self.LHV * self.efficiency
        h_out = h_in + Q_guess / max(mdot, 1.0e-9)
        return self.PR * P_in, h_out

    def guess_outlet_mdot(self, pair: tuple[str, str], mdot_in: float) -> float:
        if self.use_fuel_port:
            # Fuel arrives through its own port/node with its own mdot
            # unknown, not folded into this pair's inlet->outlet guess.
            return mdot_in
        fuel_guess = self.mdot_fuel if self.mdot_fuel is not None else MDOT_FUEL_GUESS_FRACTION * mdot_in
        return mdot_in + fuel_guess

    def residuals(self, state: "NetworkState") -> list[float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        mdot_in = state.mdot(self._inlet_node)
        mdot_fuel = self._fuel_flow(state)

        fuel_enthalpy_term = 0.0
        if self.use_fuel_port:
            _P_fuel_in, h_fuel_in = state.node(self._fuel_in_node)
            fuel_enthalpy_term = mdot_fuel * h_fuel_in

        Q = mdot_fuel * self.LHV * self.efficiency
        Q_loss = heat_loss_watts(self.heat_path, state)
        mdot_out_expected = mdot_in + mdot_fuel

        momentum_residual = P_out - self.PR * P_in
        energy_residual = mdot_out_expected * h_out - (
            mdot_in * h_in + fuel_enthalpy_term + Q - Q_loss
        )
        mass_residual = state.mdot(self._outlet_node) - mdot_out_expected
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        mdot_fuel = self._fuel_flow(state)
        return {
            "power [W]": state.mdot(self._outlet_node) * h_out - state.mdot(self._inlet_node) * h_in,
            "mdot_fuel [kg/s]": mdot_fuel,
            "PR [-]": P_out / P_in,
            "T_out [K]": state.fluid_at(self._outlet_node).temperature_ph(P_out, h_out),
            "Q_loss [W]": heat_loss_watts(self.heat_path, state),
        }
