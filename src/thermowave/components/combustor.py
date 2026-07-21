from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.components.heat_transfer import heat_loss_watts
from thermowave.core.constants import MDOT_FUEL_GUESS_FRACTION
from thermowave.fluids.base_fluid import BaseFluid
from thermowave.fluids.cantera_fluid import CanteraFluid, _CanteraCompositionFluid

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState

_REPORTED_SPECIES = ("CO2", "H2O", "O2", "N2", "CO", "NO")  # major products + common pollutants
_MOLE_FRACTION_THRESHOLD = 1.0e-6  # species below this are omitted from product_composition()


class Combustor(BaseComponent):
    """Combustor using Cantera chemical-equilibrium combustion products to
    find the outlet temperature, instead of SimpleCombustor's fixed-LHV
    heat-release model.

    Requires the optional 'cantera' extra: pip install thermowave[cantera]

    At residual-evaluation time: the inlet air (T_in, P_in from the
    network's own node state, composition `oxidizer`) is mixed by mass with
    the fuel (`fuel`, at the same T_in/P_in) in the ratio mdot_fuel/mdot_in,
    then equilibrated at constant enthalpy and pressure (`equilibrate("HP")`)
    using a Cantera mechanism (`mechanism`, default GRI-Mech 3.0) — a
    standard adiabatic-flame-temperature calculation, giving a genuinely
    composition- and dissociation-aware T_out rather than a constant-LHV
    estimate. `efficiency` (default 1.0) scales the *temperature rise* from
    that adiabatic result (T_out_actual = T_in + efficiency * (T_out_adiabatic
    - T_in)), a simplified way to represent incomplete combustion / wall
    heat loss without re-running the equilibrium at a different enthalpy.

    Composition propagation: when this Combustor's own inlet fluid (the
    network's default `fluid`, or whatever a component further upstream
    already changed it to) is itself a CanteraFluid, the reacted product
    composition is fed back into the network via outlet_fluid() — every
    component downstream (Pipe, Turbine, a heat exchanger, ...) that reads
    NetworkState.fluid_at(<its own node>) instead of `.fluid` directly then
    sees the real product mixture's density/enthalpy/cp, not plain air's.
    This is only physically consistent when the inlet fluid already shares
    Cantera's absolute (formation-enthalpy-referenced) datum — mixing that
    with e.g. IdealGasFluid's h=cp*T (referenced to h=0 at T=0K, with no
    notion of chemical/formation energy at all) would silently corrupt any
    downstream energy balance. So for any inlet fluid that isn't a
    CanteraFluid, outlet_fluid() returns None (pass-through) and this
    Combustor falls back to its original behavior: T_out is chemistry-
    informed, but the outlet enthalpy (and everything downstream) is still
    tracked through the same fluid model as upstream — "chemistry-informed
    T_out, tracked downstream as if it were still the same working fluid".
    Product *composition* is always visible either way — call
    product_composition(state) directly for the full equilibrium mole-
    fraction breakdown, or read the major species (CO2, H2O, O2, N2, CO, NO
    — whichever the mechanism produces above a trace threshold) straight
    off report_metrics() as "X_<species> [-]". SimpleCombustor is the
    fully-consistent-with-itself alternative if even chemistry-informed
    T_out is more machinery than you want.

    mdot_fuel [kg/s]: give it directly, or leave it None to drive the
    combustor by some other known quantity instead (same free-parameter /
    Setpoint-or-Controller pattern as Compressor's N — see SimpleCombustor's
    docstring for the mechanism).

    heat_path: optional Convection/Conduction/Radiation
    (thermowave.components.heat_transfer) representing heat lost through
    this combustor's liner/casing to something else, on top of whatever
    `efficiency` already models. None (the default) means no additional
    loss, unchanged from before this existed. Like Turbine/Compressor,
    this component's own energy_residual is per-unit-mass (h_out vs.
    h_out_target), so heat_path's Q(state) is added divided by the
    outlet's own expected mass flow. Since the path needs (self, "out")
    as one of its own endpoints, it can only be built after this Combustor
    already exists — pass it here if you have it, or just set
    combustor.heat_path = path afterwards.
    """

    def __init__(
        self,
        name: str,
        PR: float = 0.97,
        efficiency: float = 1.0,
        mdot_fuel: float | None = None,
        fuel: str = "CH4",
        oxidizer: str = "O2:0.21, N2:0.79",
        mechanism: str = "gri30.yaml",
        heat_path: BaseComponent | None = None,
    ):
        try:
            import cantera as ct
        except ImportError as exc:
            raise ImportError(
                "Combustor requires the 'cantera' extra: pip install thermowave[cantera]"
            ) from exc

        self._ct = ct
        self.name = name
        self.PR = PR
        self.efficiency = efficiency
        self.mdot_fuel = mdot_fuel
        self.fuel = fuel
        self.oxidizer = oxidizer
        self.mechanism = mechanism
        self.heat_path = heat_path
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"
        # ct.Solution(mechanism) parses the whole mechanism file (GRI-Mech
        # 3.0's default is 53 species/325 reactions) — loaded once here and
        # reused by every _equilibrate() call instead of re-parsing it from
        # disk on every residual evaluation. With finite-difference
        # Jacobians calling residuals() ~once per free unknown per Newton
        # iteration, re-parsing per call (the original behavior) turned a
        # single steady solve into hundreds of redundant mechanism loads.
        self._gas = ct.Solution(mechanism)
        # A second, independent Solution instance dedicated to
        # _CanteraCompositionFluid's arbitrary later (P, h)/(P, T) queries —
        # kept separate from self._gas (used only inside _equilibrate()'s
        # own air/fuel/mixture Quantity math) so a downstream component
        # reading the product fluid's properties can never interfere with
        # an equilibrium calculation in progress, or vice versa. One extra
        # one-time mechanism parse here, not one per Newton iteration.
        self._product_gas = ct.Solution(mechanism)

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def report_category(self) -> str:
        return "combustor"

    def _fuel_flow(self, state: "NetworkState") -> float:
        if self.mdot_fuel is not None:
            return self.mdot_fuel
        return state.param(f"{self.name}.mdot_fuel")

    def free_parameters(self) -> dict[str, float]:
        if self.mdot_fuel is not None:
            return {}
        return {"mdot_fuel": MDOT_FUEL_GUESS_FRACTION}

    def guess_free_parameters(
        self, fluid: "BaseFluid", P_in: float, h_in: float, mdot: float
    ) -> dict[str, float]:
        if self.mdot_fuel is not None:
            return {}
        return {"mdot_fuel": MDOT_FUEL_GUESS_FRACTION * mdot}

    def guess_outlet(self, P_in: float, h_in: float, mdot: float) -> tuple[float, float]:
        # No chemistry here — just enough of a P/h bump to keep downstream
        # free nodes' initial guess in the right order of magnitude (see
        # BaseComponent.guess_outlet's docstring).
        return self.PR * P_in, h_in

    def guess_outlet_mdot(self, pair: tuple[str, str], mdot_in: float) -> float:
        fuel_guess = self.mdot_fuel if self.mdot_fuel is not None else MDOT_FUEL_GUESS_FRACTION * mdot_in
        return mdot_in + fuel_guess

    def _equilibrate(self, T_in: float, P_in: float, mdot_air: float, mdot_fuel: float):
        ct = self._ct
        gas = self._gas

        air = ct.Quantity(gas, constant="HP")
        air.TPX = T_in, P_in, self.oxidizer
        air.mass = mdot_air

        fuel_stream = ct.Quantity(gas, constant="HP")
        fuel_stream.TPX = T_in, P_in, f"{self.fuel}:1.0"
        fuel_stream.mass = mdot_fuel

        mixture = air + fuel_stream
        mixture.equilibrate("HP")
        return mixture

    def _equilibrium_T_out(self, T_in: float, P_in: float, mdot_air: float, mdot_fuel: float) -> float:
        return float(self._equilibrate(T_in, P_in, mdot_air, mdot_fuel).T)

    def _equilibrium_mixture(self, state: "NetworkState", inlet_fluid: "BaseFluid"):
        """The Cantera equilibrium Quantity for this residual evaluation,
        memoized on `state` (see NetworkState._cache's own docstring) so
        outlet_fluid() and residuals() — both called once per Newton
        residual evaluation, see Network._resolve_node_fluid() — share one
        equilibrate() call instead of each running it independently.
        """
        key = (self.name, "mixture")
        if key not in state._cache:
            P_in, h_in = state.node(self._inlet_node)
            mdot_in = state.mdot(self._inlet_node)
            mdot_fuel = self._fuel_flow(state)
            T_in = inlet_fluid.temperature_ph(P_in, h_in)
            state._cache[key] = self._equilibrate(T_in, P_in, mdot_in, mdot_fuel)
        return state._cache[key]

    def _make_product_fluid(self, mixture) -> "_CanteraCompositionFluid":
        mass_fractions = mixture.mass_fraction_dict(threshold=1.0e-12)
        return _CanteraCompositionFluid(
            f"{self.name}.product", self._product_gas, mass_fractions, self.mechanism
        )

    def outlet_fluid(
        self, state: "NetworkState", pair: tuple[str, str], inlet_fluid: "BaseFluid"
    ) -> "BaseFluid | None":
        if pair != ("in", "out") or not isinstance(inlet_fluid, CanteraFluid):
            # Composition propagation is only physically consistent when
            # the upstream fluid already shares Cantera's absolute datum —
            # see this class's own docstring. Anything else falls back to
            # pass-through (today's original behavior).
            return None
        mixture = self._equilibrium_mixture(state, inlet_fluid)
        return self._make_product_fluid(mixture)

    def product_composition(self, state: "NetworkState") -> dict[str, float]:
        """Mole fractions of the equilibrium combustion products at the
        current (converged) inlet state and fuel flow — re-runs the same
        Cantera equilibrium calculation residuals() uses for T_out, but
        returns the full composition instead of discarding it. Only species
        with mole fraction >= _MOLE_FRACTION_THRESHOLD are included.
        """
        inlet_fluid = state.fluid_at(self._inlet_node)
        mixture = self._equilibrium_mixture(state, inlet_fluid)
        return mixture.mole_fraction_dict(threshold=_MOLE_FRACTION_THRESHOLD)

    def residuals(self, state: "NetworkState") -> list[float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        mdot_in = state.mdot(self._inlet_node)
        mdot_fuel = self._fuel_flow(state)

        inlet_fluid = state.fluid_at(self._inlet_node)
        T_in = inlet_fluid.temperature_ph(P_in, h_in)
        mixture = self._equilibrium_mixture(state, inlet_fluid)
        T_out_adiabatic = float(mixture.T)
        T_out_target = T_in + self.efficiency * (T_out_adiabatic - T_in)
        # If composition propagates (inlet_fluid is a CanteraFluid),
        # evaluate the outlet enthalpy through the same product fluid
        # outlet_fluid() hands downstream components — otherwise h_out
        # (read via inlet_fluid) and this residual would disagree about
        # what "h_out" means for the exact same (P_out, T_out) pair.
        # Falls back to inlet_fluid (pass-through) exactly when
        # outlet_fluid() does, for the same reason.
        outlet_fluid = self.outlet_fluid(state, ("in", "out"), inlet_fluid) or inlet_fluid
        h_out_target = outlet_fluid.enthalpy_pt(P_out, T_out_target)

        mdot_out_expected = mdot_in + mdot_fuel
        Q_loss = heat_loss_watts(self.heat_path, state)

        momentum_residual = P_out - self.PR * P_in
        energy_residual = h_out - h_out_target + Q_loss / mdot_out_expected
        mass_residual = state.mdot(self._outlet_node) - mdot_out_expected
        return [momentum_residual, energy_residual, mass_residual]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        mdot_fuel = self._fuel_flow(state)
        metrics = {
            "power [W]": state.mdot(self._outlet_node) * h_out - state.mdot(self._inlet_node) * h_in,
            "mdot_fuel [kg/s]": mdot_fuel,
            "PR [-]": P_out / P_in,
            "T_out [K]": state.fluid_at(self._outlet_node).temperature_ph(P_out, h_out),
            "Q_loss [W]": heat_loss_watts(self.heat_path, state),
        }
        products = self.product_composition(state)
        for species in _REPORTED_SPECIES:
            if species in products:
                metrics[f"X_{species} [-]"] = products[species]
        return metrics
