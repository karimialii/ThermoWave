from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.components.heat_transfer import heat_loss_watts
from thermowave.core.constants import MDOT_FLOOR, N_GUESS_T_FALLBACK, PA_PER_BAR
from thermowave.core.exceptions import FluidRangeError
from thermowave.maps.characteristic_map import CharacteristicMap

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid


class Compressor(BaseComponent):
    """Compressor driven by a Flownex-style (.cop) characteristic map.

    Pressure ratio and isentropic efficiency are read off the map's iso-speed
    curves (PR and eta vs non-dimensional/corrected mass flow) at the
    compressor's shaft speed N and the corrected mass flow implied by the
    current inlet state, rather than being fixed inputs as in
    SimpleCompressor. Outlet enthalpy is then found the same way as
    SimpleCompressor: isentropic rise from the ideal-gas relation, scaled up
    by 1/eta_s. gamma: give it directly, or leave it None (the default) to
    derive it from the network's own fluid model instead, via
    BaseFluid.gamma(P_in, T_in) evaluated fresh at each residual call — every
    fluid model here implements cp()/cv() (see BaseFluid.gamma()'s
    docstring), so this works for CoolProp/Cantera real-fluid models too,
    not just the constant-cp ideal-gas ones. Passing gamma directly instead
    is still useful to pin a known constant value or skip the extra
    property-model call.

    N is the shaft speed [rev/min]. Give it directly, or leave it None to
    drive the compressor by some other known quantity instead (a target
    power, a target PR, ...) — N then becomes an extra Newton unknown (via a
    free mechanical port, seeded from the map's own mid-speed) and needs a
    matching residual from somewhere else in the network to pin it down,
    e.g. a Setpoint targeting this compressor's "shaft" port, or a Shaft
    tying it to another machine's speed via
    network.connect(comp, "shaft", other, "shaft", kind="mechanical"). An
    unconnected compressor still gets its own private mechanical node (no
    different from before this port existed). This also exposes a "power"
    signal port (see provided_signal_values()) reporting the same value as
    report_metrics()["power [W]"], for a Generator/Shaft to read via
    kind="signal" instead of holding a direct object reference.

    factor_overrides: optional dict overriding any of the map file's own
    conversion factors (A_fact, B_fact, C_fact, E_fact — see
    CharacteristicMap's docstring), to calibrate the map against test data
    without editing the map file itself. Omit it (or leave a given key out)
    to use the file's own value, unchanged.

    heat_path: optional Convection/Conduction/Radiation
    (thermowave.components.heat_transfer) representing heat this
    compressor's fluid loses to something else (its own casing, ambient,
    ...). None (the default) means fully adiabatic, unchanged from before
    this existed. Since the path needs (self, "out") as one of its own
    endpoints, it can only be built after this Compressor already exists —
    pass it here if you have it, or just set it afterwards via
    comp.set(heat_path=path) (or direct attribute assignment; both work
    identically, residuals() reads the attribute either way). Q(state) > 0 (heat leaving the fluid) reduces T_out below
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
            raise ValueError(f"Compressor {name!r}: gamma must be > 1, got {gamma}")
        self.name = name
        self.map = CharacteristicMap.from_file(map_path, factor_overrides=factor_overrides)
        self.gamma = gamma
        self.N = N
        self.heat_path = heat_path
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"
        self._shaft_port = f"{name}.shaft"
        self._power_port = f"{name}.power"

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def mechanical_ports(self) -> dict[str, str]:
        return {"shaft": self._shaft_port}

    def signal_ports(self) -> dict[str, str]:
        return {"power": self._power_port}

    def shaft_sign(self) -> float:
        return -1.0  # draws power from the shaft (report_metrics()'s
        # mdot*(h_out-h_in) is positive when compressing)

    def report_category(self) -> str:
        return "turbomachinery"

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        # PR isn't known yet if N is free — a generic mid-range multiplier is
        # enough to keep downstream free nodes' initial P guess in the right
        # order of magnitude (see BaseComponent.guess_outlet's docstring).
        return 3.0 * P_in, h_in

    def fixed_mechanical_values(self) -> dict[str, float]:
        if self.N is None:
            return {}
        return {self._shaft_port: self.N}

    def free_mechanical_ports(self) -> dict[str, float]:
        if self.N is not None:
            return {}
        N_guess = self.map.mid_speed() * N_GUESS_T_FALLBACK**0.5 * 60.0
        return {self._shaft_port: N_guess}

    def guess_free_mechanical_ports(
        self, fluid: "BaseFluid", P_in: float, h_in: float, mdot: float
    ) -> dict[str, float]:
        if self.N is not None:
            return {}
        T_in = fluid.temperature_ph(P_in, h_in)
        N_guess = self.map.mid_speed() * T_in**0.5 * 60.0
        return {self._shaft_port: N_guess}

    def provided_signal_values(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._inlet_node)
        _, h_out = state.node(self._outlet_node)
        return {self._power_port: state.mdot(self._inlet_node) * (h_out - h_in)}

    def _shaft_speed(self, state: "NetworkState") -> float:
        if self.N is not None:
            return self.N
        return state.N(self._shaft_port)

    def _gamma(self, state: "NetworkState", P_in: float, T_in: float) -> float:
        if self.gamma is not None:
            return self.gamma
        return state.fluid_at(self._inlet_node).gamma(P_in, T_in)

    def _corrected_params(self, state: "NetworkState") -> tuple[float, float, float]:
        P_in, h_in = state.node(self._inlet_node)
        T_in = state.fluid_at(self._inlet_node).temperature_ph(P_in, h_in)
        mdot = state.mdot(self._inlet_node)
        if T_in <= 0.0:
            # A bad Newton iterate landing on a negative T_in would
            # otherwise make T_in**0.5 silently return a complex number
            # (Python 3's behavior for a negative float base) instead of
            # failing at the source -- surfacing as a confusing TypeError
            # far downstream instead. Same recoverable-invalid-state signal
            # a real fluid backend already raises for this trial.
            raise FluidRangeError(
                f"Compressor {self.name!r}: inlet temperature must be "
                f"positive, got T_in={T_in} K"
            )
        if mdot < 0.0:
            # Reverse/windmilling flow isn't represented in this map's
            # tabulated data -- without this check, _interpolate_1d would
            # just clamp the negative corrected flow to the map's lowest
            # tabulated point and return a normal-looking but physically
            # wrong PR/efficiency instead of flagging that this state can't
            # be evaluated.
            raise FluidRangeError(
                f"Compressor {self.name!r}: reverse flow (mdot={mdot} kg/s) "
                "is not supported by this component's characteristic map"
            )
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
        T_out_isentropic = T_in * PR ** ((gamma - 1.0) / gamma)
        h_out_isentropic = state.fluid_at(self._inlet_node).enthalpy_pt(P_out, T_out_isentropic)
        dh_actual = (h_out_isentropic - h_in) / eta_s

        mdot_in = state.mdot(self._inlet_node)
        Q_loss = heat_loss_watts(self.heat_path, state)

        momentum_residual = P_out - PR * P_in
        energy_residual = h_out - (h_in + dh_actual) + Q_loss / max(mdot_in, MDOT_FLOOR)
        mass_residual = state.mdot(self._outlet_node) - mdot_in
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        _, A, B = self._corrected_params(state)
        return {
            "power [W]": state.mdot(self._inlet_node) * (h_out - h_in),
            "eta_s [-]": self.map.efficiency(A, B),
            "PR [-]": P_out / P_in,
            "N [rev/min]": self._shaft_speed(state),
            "Q_loss [W]": heat_loss_watts(self.heat_path, state),
            "T_in [K]": state.fluid_at(self._inlet_node).temperature_ph(P_in, h_in),
            "T_out [K]": state.fluid_at(self._outlet_node).temperature_ph(P_out, h_out),
        }
