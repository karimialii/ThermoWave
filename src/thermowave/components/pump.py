from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.fluids.two_phase import require_entropy

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class Pump(BaseComponent):
    """Liquid pump: raises pressure with a small enthalpy rise, via an
    isentropic path scaled by efficiency — the liquid counterpart to a
    compressor, and the component that closes a Rankine loop's low-pressure
    side back up to boiler pressure.

    Entropy-based (works for a real liquid, and for wet steam, unlike the
    gas turbines' gamma-relation): s_in = entropy_ph(P_in, h_in);
    h_out_isentropic = enthalpy_ps(P_out, s_in); the actual enthalpy rise is
    the isentropic rise divided by eta (inefficiency means MORE work goes in
    than the ideal reversible pump). Requires a fluid exposing entropy_ph /
    enthalpy_ps (CoolPropFluid), checked at residual time.

    Specify exactly one of P_out [Pa] (absolute discharge pressure) or PR
    (P_out / P_in, > 1). 3 residuals: momentum, energy, mass.
    """

    def __init__(
        self,
        name: str,
        P_out: float | None = None,
        PR: float | None = None,
        eta: float = 0.75,
    ):
        if (P_out is None) == (PR is None):
            raise ValueError(f"Pump {name!r}: give exactly one of P_out or PR")
        if P_out is not None and P_out <= 0:
            raise ValueError(f"Pump {name!r}: P_out must be > 0, got {P_out}")
        if PR is not None and PR <= 1.0:
            raise ValueError(f"Pump {name!r}: PR must be > 1 (a pump raises pressure), got {PR}")
        if not (0.0 < eta <= 1.0):
            raise ValueError(f"Pump {name!r}: eta must be in (0, 1], got {eta}")
        self.name = name
        self.P_out = P_out
        self.PR = PR
        self.eta = eta
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def report_category(self) -> str:
        return "turbomachinery"

    def _P_out_target(self, P_in: float) -> float:
        return self.P_out if self.P_out is not None else self.PR * P_in

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        return self._P_out_target(P_in), h_in

    def residuals(self, state: "NetworkState") -> list[float]:
        fluid = state.fluid_at(self._inlet_node)
        require_entropy(fluid, f"Pump {self.name!r}")
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)

        s_in = fluid.entropy_ph(P_in, h_in)
        h_out_isentropic = fluid.enthalpy_ps(P_out, s_in)
        dh_actual = (h_out_isentropic - h_in) / self.eta

        momentum_residual = P_out - self._P_out_target(P_in)
        energy_residual = h_out - (h_in + dh_actual)
        mass_residual = state.mdot(self._outlet_node) - state.mdot(self._inlet_node)
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        return {
            "power [W]": state.mdot(self._inlet_node) * (h_out - h_in),  # work input, positive
            "eta_s [-]": self.eta,
            "PR [-]": P_out / P_in,
        }
