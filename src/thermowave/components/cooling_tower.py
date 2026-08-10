from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.fluids.humid_air import HumidAirFluid
from thermowave.fluids.psychrometrics import require_humid_air, supports_humid_air

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid


class CoolingTower(BaseComponent):
    """Evaporative cooling tower: a warm water stream is cooled by direct
    contact with an air stream, which picks up water-vapor mass (and hence
    changes its own humidity ratio W) as it does — the real, evaporative
    counterpart to modeling a cooling tower as a plain dry-air HeatExchanger
    (the SEGS benchmark's own stand-in, which has no humidity model at all).

    Four ports: water_in/water_out (liquid being cooled) and air_in/air_out
    (humid air -- air_in's fluid must be humid-air-capable, i.e. a
    HumidAirFluid; see thermowave.fluids.psychrometrics.require_humid_air()).

    Genuinely harder than Combustor's own mass-adding precedent
    (outlet_fluid() returning a new composition, mdot_out = mdot_in +
    mdot_added): Combustor's fuel mass arrives from an independent,
    unconstrained port, but here the evaporated mass must be conserved
    BETWEEN this component's own two streams -- what the water side loses,
    the air side's W must gain, exactly. One simplification falls out for
    free from HumidAirFluid's own documented mdot convention (mdot means
    DRY-AIR mass flow): dry air itself is conserved through a pure
    evaporative process, so mdot_air_out == mdot_air_in exactly (no
    guess_outlet_mdot() override needed on the air pair at all) -- only W
    changes on the air side, and only the water side needs a genuine
    mass-losing residual.

    v1 physics -- outlet air modeled as saturated (a standard, well-
    established simplified engineering idealization for a first cut; a
    Merkel/NTU-based partial-saturation effectiveness is a natural future
    extension, the same maturity path HeatExchanger itself took from fixed-
    effectiveness to UA/NTU -- see target_RH_out below for the non-breaking
    hook a future version would use):

      range [K]: water_in's temperature minus water_out's target
          (mirrors Condenser's direct-target-state style, not an
          effectiveness calc) -- h_water_out_target = water_fluid.enthalpy_pt(
          P_water_out, T_water_in - range).
      W_air_out [kg/kg]: a genuinely free Newton unknown (like Combustor's
          own mdot_fuel), closed by a saturation residual pinning the
          outlet air's own relative humidity to target_RH_out (default 1.0,
          fully saturated -- the idealization above). outlet_fluid()
          returns a fresh HumidAirFluid built from this resolved W_air_out
          each residual evaluation, exactly the Combustor pattern for
          propagating a resolved, per-iteration composition downstream.
      Overall control-volume energy balance (not a two-step "duty then
          push the same duty across" the way Condenser/FeedwaterHeater
          are -- there's no single well-defined "Q" here once evaporated
          mass carries its own enthalpy across the boundary with it, and
          HumidAirFluid's own H is already, by construction, a dry-air-
          basis enthalpy that correctly includes whatever water-vapor
          content each humid-air state carries):
          mdot_water_in*h_water_in + mdot_dry_air*h_air_in ==
          mdot_water_out*h_water_out_target + mdot_dry_air*h_air_out

    7 residuals: water momentum/energy/mass, air momentum/energy/mass, and
    the saturation residual closing W_air_out.
    """

    def __init__(
        self,
        name: str,
        range: float,
        PR_water: float = 1.0,
        PR_air: float = 1.0,
        target_RH_out: float = 1.0,
        W_air_out_guess: float = 0.02,
    ):
        if range <= 0.0:
            raise ValueError(f"CoolingTower {name!r}: range must be > 0, got {range}")
        if not (0.0 < PR_water <= 1.0):
            raise ValueError(f"CoolingTower {name!r}: PR_water must be in (0, 1], got {PR_water}")
        if not (0.0 < PR_air <= 1.0):
            raise ValueError(f"CoolingTower {name!r}: PR_air must be in (0, 1], got {PR_air}")
        if not (0.0 < target_RH_out <= 1.0):
            raise ValueError(
                f"CoolingTower {name!r}: target_RH_out must be in (0, 1], got {target_RH_out}"
            )
        self.name = name
        self.range = range
        self.PR_water = PR_water
        self.PR_air = PR_air
        self.target_RH_out = target_RH_out
        self.W_air_out_guess = W_air_out_guess
        self._water_in_node = f"{name}.water_in"
        self._water_out_node = f"{name}.water_out"
        self._air_in_node = f"{name}.air_in"
        self._air_out_node = f"{name}.air_out"

    def ports(self) -> dict[str, str]:
        return {
            "water_in": self._water_in_node,
            "water_out": self._water_out_node,
            "air_in": self._air_in_node,
            "air_out": self._air_out_node,
        }

    def warm_start_pairs(self) -> list[tuple[str, str]]:
        return [("water_in", "water_out"), ("air_in", "air_out")]

    def guess_outlet_for_pair(
        self, pair: tuple[str, str], P_in: float, h_in: float, mdot: float
    ) -> tuple[float, float]:
        if pair == ("water_in", "water_out"):
            return self.PR_water * P_in, h_in - 2.0e5  # a plausible sensible drop
        if pair == ("air_in", "air_out"):
            return self.PR_air * P_in, h_in + 6.5e4  # a plausible sensible+latent rise
        return P_in, h_in

    def free_parameters(self) -> dict[str, float]:
        return {"W_air_out": self.W_air_out_guess}

    def _air_out_fluid(self, state: "NetworkState") -> HumidAirFluid:
        W_air_out = state.param(f"{self.name}.W_air_out")
        return HumidAirFluid(name=f"{self.name}.air_out", W=W_air_out)

    def outlet_fluid(
        self, state: "NetworkState", pair: tuple[str, str], inlet_fluid: "BaseFluid"
    ) -> "BaseFluid | None":
        if pair != ("air_in", "air_out"):
            return None
        return self._air_out_fluid(state)

    def residuals(self, state: "NetworkState") -> list[float]:
        air_fluid = state.fluid_at(self._air_in_node)
        require_humid_air(air_fluid, f"CoolingTower {self.name!r}")
        assert supports_humid_air(air_fluid)
        water_fluid = state.fluid_at(self._water_in_node)

        P_water_in, h_water_in = state.node(self._water_in_node)
        P_water_out, h_water_out = state.node(self._water_out_node)
        P_air_in, h_air_in = state.node(self._air_in_node)
        P_air_out, h_air_out = state.node(self._air_out_node)

        mdot_water_in = state.mdot(self._water_in_node)
        mdot_dry_air = state.mdot(self._air_in_node)

        T_water_in = water_fluid.temperature_ph(P_water_in, h_water_in)
        h_water_out_target = water_fluid.enthalpy_pt(P_water_out, T_water_in - self.range)

        air_out_fluid = self._air_out_fluid(state)
        T_air_out = air_out_fluid.temperature_ph(P_air_out, h_air_out)

        mdot_evap = mdot_dry_air * (air_out_fluid.W - air_fluid.W)
        mdot_water_out_target = mdot_water_in - mdot_evap

        water_momentum = P_water_out - self.PR_water * P_water_in
        water_energy = h_water_out - h_water_out_target
        water_mass = state.mdot(self._water_out_node) - mdot_water_out_target

        air_momentum = P_air_out - self.PR_air * P_air_in
        # Overall control-volume energy balance -- see this class's own
        # docstring for why this is one shared balance rather than a
        # Condenser-style "compute Q on one side, push it across."
        air_energy = (
            mdot_water_in * h_water_in
            + mdot_dry_air * h_air_in
            - mdot_water_out_target * h_water_out_target
            - mdot_dry_air * h_air_out
        )
        air_mass = state.mdot(self._air_out_node) - mdot_dry_air

        # W_target(T_air_out) via RH as an INPUT, not relative_humidity_pt()
        # (RH as an output) -- see HumidAirFluid.humidity_ratio_at_rh's own
        # docstring: a Newton trial iterate legitimately passes through
        # marginally-supersaturated states on the way to convergence, which
        # would make relative_humidity_pt() hard-fail outright.
        W_target = air_out_fluid.humidity_ratio_at_rh(P_air_out, T_air_out, self.target_RH_out)
        saturation_residual = air_out_fluid.W - W_target

        return [
            water_momentum, water_energy, water_mass,
            air_momentum, air_energy, air_mass,
            saturation_residual,
        ]

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        air_fluid = state.fluid_at(self._air_in_node)
        require_humid_air(air_fluid, f"CoolingTower {self.name!r}")
        assert supports_humid_air(air_fluid)
        water_fluid = state.fluid_at(self._water_in_node)

        P_water_in, h_water_in = state.node(self._water_in_node)
        P_water_out, h_water_out = state.node(self._water_out_node)
        P_air_out, h_air_out = state.node(self._air_out_node)

        T_water_in = water_fluid.temperature_ph(P_water_in, h_water_in)
        T_water_out = water_fluid.temperature_ph(P_water_out, h_water_out)

        air_out_fluid = self._air_out_fluid(state)
        T_air_out = air_out_fluid.temperature_ph(P_air_out, h_air_out)

        mdot_dry_air = state.mdot(self._air_in_node)
        mdot_evap = mdot_dry_air * (air_out_fluid.W - air_fluid.W)

        return {
            "range [K]": T_water_in - T_water_out,
            "T_water_in [K]": T_water_in,
            "T_water_out [K]": T_water_out,
            "T_air_out [K]": T_air_out,
            "W_air_in [-]": air_fluid.W,
            "W_air_out [-]": air_out_fluid.W,
            "RH_air_out [-]": air_out_fluid.relative_humidity_pt(P_air_out, T_air_out),
            "mdot_evap [kg/s]": mdot_evap,
        }

    # No report_category() override: no existing categorized table's column
    # layout fits a mass-transfer component with both a "range" and a
    # humidity ratio to report -- report_metrics() above still works
    # standalone, just outside print_report()'s categorized-tables section.
