from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.components.heat_transfer import heat_loss_watts
from thermowave.core.constants import N_GUESS_T_FALLBACK, PA_PER_BAR
from thermowave.maps.characteristic_map import CharacteristicMap

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid


class Turbine(BaseComponent):
    """Turbine driven by a Flownex-style (.tur) characteristic map.

    Same file format and corrected groups as Compressor (A = corrected
    speed, B = corrected mass flow, both referenced to turbine *inlet*
    conditions, matching the map's own P03/T03 = inlet convention). The map's
    "pressure ratio" column C is the turbine's expansion ratio P_in/P_out
    (>1, per the file's own "C = P03/P04" comment — inlet over outlet, the
    opposite sense from the compressor map, where the same column reads as
    P_out/P_in).

    PR and isentropic efficiency are looked up at the turbine's rotational
    speed N and the corrected mass flow implied by the current inlet state,
    rather than being fixed inputs as in SimpleTurbine. Outlet enthalpy is
    then found the same way as SimpleTurbine: isentropic drop from the
    ideal-gas relation, scaled down by eta_s. gamma: give it directly, or
    leave it None (the default) to derive it from the network's own fluid
    model instead, via BaseFluid.gamma(P_in, T_in) evaluated fresh at each
    residual call — every fluid model here implements cp()/cv() (see
    BaseFluid.gamma()'s docstring), so this works for CoolProp/Cantera
    real-fluid models too, not just the constant-cp ideal-gas ones. Passing
    gamma directly instead is still useful to pin a known constant value or
    skip the extra property-model call.

    N is the shaft speed [rev/min]. Give it directly, or leave it None to
    drive the turbine by some other known quantity instead (a target power,
    a target PR, ...) — N then becomes an extra Newton unknown (via
    free_parameters(), seeded from the map's own mid-speed) and needs a
    matching residual from somewhere else in the network to pin it down,
    e.g. a Setpoint component tying report_metrics()["power [W]"] to a
    target. There's still no shaft/mechanical network connection: N is
    either given directly or solved purely from residuals contributed
    elsewhere, not coupled to any other component's speed.

    factor_overrides: optional dict overriding any of the map file's own
    conversion factors (A_fact, B_fact, C_fact, E_fact — see
    CharacteristicMap's docstring), to calibrate the map against test data
    without editing the map file itself. Omit it (or leave a given key out)
    to use the file's own value, unchanged.

    heat_path: optional Convection/Conduction/Radiation
    (thermowave.components.heat_transfer) representing heat this turbine's
    fluid loses to something else (its own casing, ambient, ...). None
    (the default) means fully adiabatic, unchanged from before this
    existed. Since the path needs (self, "out") as one of its own
    endpoints, it can only be built after this Turbine already exists —
    pass it here if you have it, or just set turb.heat_path = path
    afterwards; both work identically, residuals() reads the attribute
    either way. Q(state) > 0 (heat leaving the fluid) reduces T_out below
    what the map/efficiency alone would give, same sign convention as
    Pipe's own heat_loss.
    """

    def __init__(
        self,
        name: str,
        map_path: str,
        gamma: float | None = None,
        N: float | None = None,
        factor_overrides: dict[str, float] | None = None,
        heat_path: BaseComponent | None = None,
    ):
        if gamma is not None and gamma <= 1.0:
            raise ValueError(f"Turbine {name!r}: gamma must be > 1, got {gamma}")
        self.name = name
        self.map = CharacteristicMap.from_file(map_path, factor_overrides=factor_overrides)
        self.gamma = gamma
        self.N = N
        self.heat_path = heat_path
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def report_category(self) -> str:
        return "turbomachinery"

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        # PR isn't known yet if N is free — a generic mid-range multiplier is
        # enough to keep downstream free nodes' initial P guess in the right
        # order of magnitude (see BaseComponent.guess_outlet's docstring).
        return P_in / 2.0, h_in

    def free_parameters(self) -> dict[str, float]:
        if self.N is not None:
            return {}
        N_guess = self.map.mid_speed() * N_GUESS_T_FALLBACK**0.5 * 60.0
        return {"N": N_guess}

    def guess_free_parameters(
        self, fluid: "BaseFluid", P_in: float, h_in: float, mdot: float
    ) -> dict[str, float]:
        if self.N is not None:
            return {}
        T_in = fluid.temperature_ph(P_in, h_in)
        N_guess = self.map.mid_speed() * T_in**0.5 * 60.0
        return {"N": N_guess}

    def _shaft_speed(self, state: "NetworkState") -> float:
        if self.N is not None:
            return self.N
        return state.param(f"{self.name}.N")

    def _gamma(self, state: "NetworkState", P_in: float, T_in: float) -> float:
        if self.gamma is not None:
            return self.gamma
        return state.fluid_at(self._inlet_node).gamma(P_in, T_in)

    def _corrected_params(self, state: "NetworkState") -> tuple[float, float, float]:
        P_in, h_in = state.node(self._inlet_node)
        T_in = state.fluid_at(self._inlet_node).temperature_ph(P_in, h_in)
        mdot = state.mdot(self._inlet_node)
        N = self._shaft_speed(state)
        A = (N / 60.0) / T_in**0.5
        B = mdot * T_in**0.5 / (P_in / PA_PER_BAR)
        return T_in, A, B

    def residuals(self, state: "NetworkState") -> list[float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)

        T_in, A, B = self._corrected_params(state)
        PR = self.map.pressure_ratio(A, B)
        eta_s = self.map.efficiency(A, B)

        gamma = self._gamma(state, P_in, T_in)
        T_out_isentropic = T_in * (1.0 / PR) ** ((gamma - 1.0) / gamma)
        h_out_isentropic = state.fluid_at(self._inlet_node).enthalpy_pt(P_out, T_out_isentropic)
        dh_actual = eta_s * (h_in - h_out_isentropic)

        mdot_in = state.mdot(self._inlet_node)
        Q_loss = heat_loss_watts(self.heat_path, state)

        momentum_residual = P_in - PR * P_out
        energy_residual = h_out - (h_in - dh_actual) + Q_loss / mdot_in
        mass_residual = state.mdot(self._outlet_node) - mdot_in
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        _, A, B = self._corrected_params(state)
        return {
            "power [W]": state.mdot(self._inlet_node) * (h_in - h_out),
            "eta_s [-]": self.map.efficiency(A, B),
            "PR [-]": P_in / P_out,
            "N [rev/min]": self._shaft_speed(state),
            "Q_loss [W]": heat_loss_watts(self.heat_path, state),
        }
