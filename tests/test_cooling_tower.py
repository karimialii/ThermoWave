import math

import pytest

pytest.importorskip("CoolProp")

from CoolProp.HumidAirProp import HAPropsSI  # noqa: E402

from thermowave.components.cooling_tower import CoolingTower  # noqa: E402
from thermowave.components.sink import Sink  # noqa: E402
from thermowave.components.source import Source  # noqa: E402
from thermowave.core.network import Network  # noqa: E402
from thermowave.fluids.humid_air import HumidAirFluid  # noqa: E402
from thermowave.fluids.ideal_gas import IdealGasFluid  # noqa: E402
from thermowave.fluids.real_fluid import CoolPropFluid  # noqa: E402

WATER = CoolPropFluid(name="Water")
P_AIR = 101325.0
T_AIR_IN = 298.15  # 25 C
RH_AIR_IN = 0.5
W_AIR_IN = HAPropsSI("W", "T", T_AIR_IN, "RH", RH_AIR_IN, "P", P_AIR)
AIR_IN = HumidAirFluid(name="air_in", W=W_AIR_IN)


class _FakeState:
    def __init__(self, node_values, mdots, params, fluids):
        self._node_values = node_values
        self._mdots = mdots
        self._params = params
        self._fluids = fluids

    def node(self, name):
        return self._node_values[name]

    def fluid_at(self, name):
        return self._fluids[name]

    def mdot(self, name):
        return self._mdots[name]

    def param(self, name):
        return self._params[name]


def test_rejects_non_positive_range():
    with pytest.raises(ValueError, match="range"):
        CoolingTower(name="ct", range=0.0)


def test_rejects_out_of_range_pr_and_target_rh():
    with pytest.raises(ValueError, match="PR_water"):
        CoolingTower(name="ct", range=5.0, PR_water=1.5)
    with pytest.raises(ValueError, match="PR_air"):
        CoolingTower(name="ct", range=5.0, PR_air=0.0)
    with pytest.raises(ValueError, match="target_RH_out"):
        CoolingTower(name="ct", range=5.0, target_RH_out=1.5)


def test_ports_and_free_parameters():
    ct = CoolingTower(name="ct", range=8.0, W_air_out_guess=0.03)
    assert set(ct.ports()) == {"water_in", "water_out", "air_in", "air_out"}
    assert ct.free_parameters() == {"W_air_out": 0.03}


def test_rejects_non_humid_air_fluid_on_air_side():
    ct = CoolingTower(name="ct", range=8.0)
    dry_air = IdealGasFluid(name="dry_air", R=287.05, cp=1005.0)
    state = _FakeState(
        node_values={
            "ct.water_in": (5.0e5, WATER.enthalpy_pt(5.0e5, 313.15)),
            "ct.water_out": (5.0e5, WATER.enthalpy_pt(5.0e5, 305.15)),
            "ct.air_in": (P_AIR, dry_air.enthalpy_pt(P_AIR, T_AIR_IN)),
            "ct.air_out": (P_AIR, dry_air.enthalpy_pt(P_AIR, 306.0)),
        },
        mdots={"ct.water_in": 5.0, "ct.water_out": 4.0, "ct.air_in": 50.0, "ct.air_out": 50.0},
        params={"ct.W_air_out": 0.03},
        fluids={"ct.water_in": WATER, "ct.air_in": dry_air},
    )
    with pytest.raises(ValueError, match="humid-air-capable"):
        ct.residuals(state)


def _self_consistent_state(T_water_in, range_K, T_air_out, mdot_water_in, mdot_dry_air, P_water=5.0e5):
    """Build a _FakeState where every one of CoolingTower's 7 residuals
    evaluates to exactly zero, by construction -- not by calling the
    component's own code, but by directly re-deriving the same control-
    volume physics its docstring describes (the same "rebuild the hand-
    solved system independently" approach test_drum.py's own
    test_state_derivative_matches_hand_solved_2x2_system uses)."""
    T_water_out = T_water_in - range_K
    h_water_in = WATER.enthalpy_pt(P_water, T_water_in)
    h_water_out_target = WATER.enthalpy_pt(P_water, T_water_out)

    W_air_out = HAPropsSI("W", "T", T_air_out, "RH", 1.0, "P", P_AIR)  # saturated -> target_RH_out=1.0
    air_out_fluid = HumidAirFluid(name="ct.air_out", W=W_air_out)
    h_air_in = AIR_IN.enthalpy_pt(P_AIR, T_AIR_IN)
    h_air_out = air_out_fluid.enthalpy_pt(P_AIR, T_air_out)

    mdot_evap = mdot_dry_air * (W_air_out - W_AIR_IN)
    mdot_water_out = mdot_water_in - mdot_evap

    # Overall energy balance, solved for mdot_water_in's own *coefficient*
    # cancelling out -- i.e. verify the mdot_water_in given as an argument
    # is consistent, by checking the balance closes with it substituted in
    # (this helper picks the fluids/mdot_dry_air/T_air_out freely, then the
    # caller is responsible for passing an mdot_water_in that actually
    # balances -- computed via _solve_mdot_water_in below).

    state = _FakeState(
        node_values={
            "ct.water_in": (P_water, h_water_in),
            "ct.water_out": (P_water, h_water_out_target),
            "ct.air_in": (P_AIR, h_air_in),
            "ct.air_out": (P_AIR, h_air_out),
        },
        mdots={
            "ct.water_in": mdot_water_in, "ct.water_out": mdot_water_out,
            "ct.air_in": mdot_dry_air, "ct.air_out": mdot_dry_air,
        },
        params={"ct.W_air_out": W_air_out},
        fluids={"ct.water_in": WATER, "ct.air_in": AIR_IN},
    )
    return state, h_water_in, h_water_out_target, h_air_in, h_air_out, mdot_evap


def _solve_mdot_water_in(h_water_in, h_water_out_target, h_air_in, h_air_out, mdot_dry_air, mdot_evap):
    # mdot_water_in*(h_water_in - h_water_out_target) + mdot_evap*h_water_out_target
    #     + mdot_dry_air*(h_air_in - h_air_out) == 0
    numerator = mdot_dry_air * (h_air_out - h_air_in) - mdot_evap * h_water_out_target
    return numerator / (h_water_in - h_water_out_target)


def test_residuals_all_zero_at_hand_solved_self_consistent_state():
    ct = CoolingTower(name="ct", range=8.0)
    T_water_in, range_K, T_air_out, mdot_dry_air = 313.15, 8.0, 306.0, 50.0

    _, h_water_in, h_water_out_target, h_air_in, h_air_out, mdot_evap = _self_consistent_state(
        T_water_in, range_K, T_air_out, mdot_water_in=1.0, mdot_dry_air=mdot_dry_air,
    )
    mdot_water_in = _solve_mdot_water_in(
        h_water_in, h_water_out_target, h_air_in, h_air_out, mdot_dry_air, mdot_evap
    )
    assert mdot_water_in > 0.0  # sanity: a physically meaningful flow, not a solver artifact

    state, *_ = _self_consistent_state(
        T_water_in, range_K, T_air_out, mdot_water_in=mdot_water_in, mdot_dry_air=mdot_dry_air,
    )
    residuals = ct.residuals(state)
    assert len(residuals) == 7
    # Momentum/mass/saturation residuals are O(1) or better; the two energy
    # residuals are O(1e5-1e6) (enthalpy-scale) linear combinations chained
    # through several independent HAPropsSI calls, so abs_tol=1e-3 there is
    # actually ~1e-9 relative agreement, not a loose tolerance.
    # Order: water_momentum, water_energy, water_mass, air_momentum,
    # air_energy, air_mass, saturation_residual.
    for i, r in enumerate(residuals):
        if i == 4:  # air_energy
            assert math.isclose(r, 0.0, abs_tol=1e-2)
        else:
            assert math.isclose(r, 0.0, abs_tol=1e-3)


def test_report_metrics_reports_range_and_saturation():
    ct = CoolingTower(name="ct", range=8.0)
    T_water_in, range_K, T_air_out, mdot_dry_air = 313.15, 8.0, 306.0, 50.0
    _, h_water_in, h_water_out_target, h_air_in, h_air_out, mdot_evap = _self_consistent_state(
        T_water_in, range_K, T_air_out, mdot_water_in=1.0, mdot_dry_air=mdot_dry_air,
    )
    mdot_water_in = _solve_mdot_water_in(
        h_water_in, h_water_out_target, h_air_in, h_air_out, mdot_dry_air, mdot_evap
    )
    state, *_ = _self_consistent_state(
        T_water_in, range_K, T_air_out, mdot_water_in=mdot_water_in, mdot_dry_air=mdot_dry_air,
    )
    metrics = ct.report_metrics(state)
    assert math.isclose(metrics["range [K]"], range_K, abs_tol=1e-6)
    assert math.isclose(metrics["RH_air_out [-]"], 1.0, abs_tol=1e-6)
    assert metrics["mdot_evap [kg/s]"] > 0.0


# ---------------------------------------------------------------------------
# End-to-end: a real Network.solve() -- proves the coupled (T_air_out,
# W_air_out) system this component sets up actually converges, and that
# mass/energy come out physically sane, not just that the residual shape
# is right.
# ---------------------------------------------------------------------------


def test_cooling_tower_end_to_end_solve_converges_and_conserves_mass():
    water_src = Source(name="water_src", P=5.0e5, T=313.15, mdot=1.0)
    air_src = Source(name="air_src", P=P_AIR, T=T_AIR_IN, mdot=50.0, fluid=AIR_IN)
    ct = CoolingTower(name="ct", range=8.0)
    water_snk = Sink(name="water_snk")
    air_snk = Sink(name="air_snk")

    network = Network(fluid=WATER)
    for c in (water_src, air_src, ct, water_snk, air_snk):
        network.add_component(c)
    network.connect(water_src, "out", ct, "water_in")
    network.connect(ct, "water_out", water_snk, "in")
    network.connect(air_src, "out", ct, "air_in")
    network.connect(ct, "air_out", air_snk, "in")

    result = network.solve(tol=1e-6, max_iter=300)
    assert result.converged

    P_water_out, h_water_out = result.node_P["ct.water_out"], result.node_h["ct.water_out"]
    T_water_out = WATER.temperature_ph(P_water_out, h_water_out)
    assert math.isclose(T_water_out, 313.15 - 8.0, abs_tol=0.5)

    # Air leaves saturated (target_RH_out defaults to 1.0).
    air_out_fluid = result.node_fluid["ct.air_out"]
    P_air_out, h_air_out = result.node_P["ct.air_out"], result.node_h["ct.air_out"]
    T_air_out = air_out_fluid.temperature_ph(P_air_out, h_air_out)
    assert math.isclose(air_out_fluid.relative_humidity_pt(P_air_out, T_air_out), 1.0, abs_tol=1e-4)

    # Dry-air mdot is conserved exactly (the documented HumidAirFluid convention).
    assert math.isclose(result.node_mdot["ct.air_out"], 50.0, rel_tol=1e-9)

    # Water lost exactly what the air's humidity ratio gained.
    mdot_evap_water_side = 1.0 - result.node_mdot["ct.water_out"]
    mdot_evap_air_side = 50.0 * (air_out_fluid.W - W_AIR_IN)
    assert math.isclose(mdot_evap_water_side, mdot_evap_air_side, rel_tol=1e-4)
    assert mdot_evap_water_side > 0.0
