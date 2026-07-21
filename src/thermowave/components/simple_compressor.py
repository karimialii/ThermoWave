from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.components.heat_transfer import heat_loss_watts

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class SimpleCompressor(BaseComponent):
    """Analytic compressor: fixed pressure ratio + isentropic efficiency.

    No performance map — the ideal (isentropic) outlet temperature is the
    textbook ideal-gas relation T2s = T1 * PR**((gamma-1)/gamma); the actual
    enthalpy rise is that isentropic rise divided by eta_s (inefficiency
    shows up as extra enthalpy rise, i.e. extra heating of the fluid, the
    standard adiabatic-irreversible-compression model). gamma: give it
    directly, or leave it None (the default) to derive it from the
    network's own fluid model instead, via BaseFluid.gamma(P_in, T_in)
    evaluated fresh at each residual call — every fluid model here
    implements cp()/cv() (see BaseFluid.gamma()'s docstring), so this works
    for CoolProp/Cantera real-fluid models too, not just the constant-cp
    ideal-gas ones. Passing gamma directly instead is still useful to pin a
    known constant value or skip the extra property-model call. A map-based
    Compressor (Flownex-style iso-speed characteristic curves) is a
    separate, more detailed component.

    heat_path: optional Convection/Conduction/Radiation
    (thermowave.components.heat_transfer) representing heat this
    compressor's fluid loses to something else. None (the default) means
    fully adiabatic, unchanged from before this existed — see Compressor's
    docstring for the full explanation (same mechanism, same sign
    convention as Pipe's heat_loss).
    """

    def __init__(
        self,
        name: str,
        PR: float,
        eta_s: float,
        gamma: float | None = None,
        heat_path: BaseComponent | None = None,
    ):
        if PR <= 0:
            raise ValueError(f"SimpleCompressor {name!r}: PR must be > 0, got {PR}")
        if not (0.0 < eta_s <= 1.0):
            raise ValueError(f"SimpleCompressor {name!r}: eta_s must be in (0, 1], got {eta_s}")
        if gamma is not None and gamma <= 1.0:
            raise ValueError(f"SimpleCompressor {name!r}: gamma must be > 1, got {gamma}")
        self.name = name
        self.PR = PR
        self.eta_s = eta_s
        self.gamma = gamma
        self.heat_path = heat_path
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def report_category(self) -> str:
        return "turbomachinery"

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        return self.PR * P_in, h_in

    def residuals(self, state: "NetworkState") -> list[float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)

        inlet_fluid = state.fluid_at(self._inlet_node)
        T_in = inlet_fluid.temperature_ph(P_in, h_in)
        gamma = self.gamma if self.gamma is not None else inlet_fluid.gamma(P_in, T_in)
        T_out_isentropic = T_in * self.PR ** ((gamma - 1.0) / gamma)
        h_out_isentropic = inlet_fluid.enthalpy_pt(P_out, T_out_isentropic)
        dh_actual = (h_out_isentropic - h_in) / self.eta_s

        mdot_in = state.mdot(self._inlet_node)
        Q_loss = heat_loss_watts(self.heat_path, state)

        momentum_residual = P_out - self.PR * P_in
        energy_residual = h_out - (h_in + dh_actual) + Q_loss / mdot_in
        mass_residual = state.mdot(self._outlet_node) - mdot_in
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        return {
            "power [W]": state.mdot(self._inlet_node) * (h_out - h_in),
            "eta_s [-]": self.eta_s,
            "PR [-]": self.PR,
            "Q_loss [W]": heat_loss_watts(self.heat_path, state),
        }
