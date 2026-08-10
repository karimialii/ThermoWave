import math

import pytest

pytest.importorskip("CoolProp")

from CoolProp.HumidAirProp import HAPropsSI  # noqa: E402

from thermowave.components.pipe import Pipe  # noqa: E402
from thermowave.components.sink import Sink  # noqa: E402
from thermowave.components.source import Source  # noqa: E402
from thermowave.core.network import Network  # noqa: E402
from thermowave.fluids.humid_air import R_DA, R_V, HumidAirFluid  # noqa: E402
from thermowave.fluids.ideal_gas import IdealGasFluid  # noqa: E402
from thermowave.fluids.psychrometrics import require_humid_air, supports_humid_air  # noqa: E402

# 20 C, 50% RH, 1 atm -- the standard textbook reference point.
P_REF = 101325.0
T_REF = 293.15
W_REF = HAPropsSI("W", "T", T_REF, "RH", 0.5, "P", P_REF)


# ---------------------------------------------------------------------------
# Constructor validation
# ---------------------------------------------------------------------------


def test_rejects_negative_w():
    with pytest.raises(ValueError, match="W must be >= 0"):
        HumidAirFluid(name="moist_air", W=-0.001)


def test_rejects_absurdly_high_w():
    with pytest.raises(ValueError, match="not physically reachable"):
        HumidAirFluid(name="moist_air", W=5.0)


def test_from_relative_humidity_derives_expected_w():
    fluid = HumidAirFluid.from_relative_humidity(name="moist_air", T=T_REF, RH=0.5, P=P_REF)
    assert math.isclose(fluid.W, W_REF, rel_tol=1e-9)
    # Cross-checked against the well-known ~0.0073 kg/kg textbook value for
    # 20 C / 50% RH / 1 atm, not just internal self-consistency.
    assert math.isclose(fluid.W, 0.00729, rel_tol=2e-2)


# ---------------------------------------------------------------------------
# The 5 BaseFluid methods, against direct HAPropsSI calls as ground truth
# ---------------------------------------------------------------------------


def test_enthalpy_pt_matches_direct_hapropssi_call():
    fluid = HumidAirFluid(name="moist_air", W=W_REF)
    h = fluid.enthalpy_pt(P_REF, T_REF)
    expected = HAPropsSI("H", "T", T_REF, "P", P_REF, "W", W_REF)
    assert math.isclose(h, expected, rel_tol=1e-9)


def test_temperature_ph_matches_direct_hapropssi_call():
    fluid = HumidAirFluid(name="moist_air", W=W_REF)
    h = HAPropsSI("H", "T", T_REF, "P", P_REF, "W", W_REF)
    T = fluid.temperature_ph(P_REF, h)
    expected = HAPropsSI("T", "P", P_REF, "H", h, "W", W_REF)
    assert math.isclose(T, expected, rel_tol=1e-9)
    assert math.isclose(T, T_REF, rel_tol=1e-6)


def test_temperature_ph_round_trips_enthalpy_pt():
    fluid = HumidAirFluid(name="moist_air", W=W_REF)
    for T in (280.0, 293.15, 310.0, 330.0):
        h = fluid.enthalpy_pt(P_REF, T)
        T_back = fluid.temperature_ph(P_REF, h)
        assert math.isclose(T_back, T, rel_tol=1e-9)


def test_density_ph_matches_inverse_of_vha():
    fluid = HumidAirFluid(name="moist_air", W=W_REF)
    h = fluid.enthalpy_pt(P_REF, T_REF)
    rho = fluid.density_ph(P_REF, h)
    Vha = HAPropsSI("Vha", "T", T_REF, "P", P_REF, "W", W_REF)
    assert math.isclose(rho, 1.0 / Vha, rel_tol=1e-9)


def test_cp_matches_cp_ha_reference():
    fluid = HumidAirFluid(name="moist_air", W=W_REF)
    cp = fluid.cp(P_REF, T_REF)
    expected = HAPropsSI("cp_ha", "T", T_REF, "P", P_REF, "W", W_REF)
    assert math.isclose(cp, expected, rel_tol=1e-9)
    # cp_ha (per kg dry air) and cp (per kg humid air) are genuinely
    # different keys -- confirm this implementation picked the dry-air one,
    # not the humid-air one, matching the documented mdot convention.
    cp_humid_air_basis = HAPropsSI("cp", "T", T_REF, "P", P_REF, "W", W_REF)
    assert not math.isclose(cp, cp_humid_air_basis, rel_tol=1e-6)


def test_cv_is_less_than_cp_and_documents_approximation():
    fluid = HumidAirFluid(name="moist_air", W=W_REF)
    cp = fluid.cp(P_REF, T_REF)
    cv = fluid.cv(P_REF, T_REF)
    assert cv < cp
    assert math.isclose(cp - cv, R_DA + W_REF * R_V, rel_tol=1e-9)


def test_gamma_default_uses_cp_over_cv():
    fluid = HumidAirFluid(name="moist_air", W=W_REF)
    assert math.isclose(
        fluid.gamma(P_REF, T_REF), fluid.cp(P_REF, T_REF) / fluid.cv(P_REF, T_REF), rel_tol=1e-12
    )


# ---------------------------------------------------------------------------
# supports_humid_air() / require_humid_air()
# ---------------------------------------------------------------------------


def test_supports_humid_air_true_for_humid_air_fluid():
    assert supports_humid_air(HumidAirFluid(name="moist_air", W=W_REF))


def test_supports_humid_air_false_for_ideal_gas_fluid():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    assert not supports_humid_air(air)


def test_require_humid_air_raises_helpful_error_for_unsupported_fluid():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    with pytest.raises(ValueError, match="needs a humid-air-capable fluid"):
        require_humid_air(air, "SomeComponent")


# ---------------------------------------------------------------------------
# Psychrometric extras: relative_humidity_pt / wet_bulb_pt / dew_point_pt
# ---------------------------------------------------------------------------


def test_relative_humidity_pt_recovers_construction_rh():
    fluid = HumidAirFluid.from_relative_humidity(name="moist_air", T=T_REF, RH=0.5, P=P_REF)
    rh = fluid.relative_humidity_pt(P_REF, T_REF)
    assert math.isclose(rh, 0.5, rel_tol=1e-6)


def test_wet_bulb_and_dew_point_bracket_dry_bulb_for_unsaturated_air():
    fluid = HumidAirFluid(name="moist_air", W=W_REF)
    T_wb = fluid.wet_bulb_pt(P_REF, T_REF)
    T_dp = fluid.dew_point_pt(P_REF, T_REF)
    # Physical sanity for unsaturated air: dew point <= wet bulb <= dry bulb,
    # strict since W_REF (50% RH) is genuinely unsaturated.
    assert T_dp < T_wb < T_REF


# ---------------------------------------------------------------------------
# End-to-end: Source(fluid=HumidAirFluid(...)) through a real Network.solve(),
# proving this needs zero network.py changes -- the composition-propagation
# machinery (fluid_seed()/_can_change_composition()) already handles a
# second fluid, verified in tests/test_fluid_propagation.py for a generic
# IdealGasFluid swap; this proves the same path for HumidAirFluid.
# ---------------------------------------------------------------------------


def test_source_humid_air_through_pipe_sensible_flow_end_to_end():
    fluid = HumidAirFluid(name="moist_air", W=W_REF)
    src = Source(name="src", P=P_REF, T=T_REF, mdot=1.0, fluid=fluid)
    pipe = Pipe(name="pipe", L=1.0, D=0.1, f=0.02)
    snk = Sink(name="snk")

    network = Network(fluid=IdealGasFluid(name="dry_air", R=287.05, cp=1005.0))
    for c in (src, pipe, snk):
        network.add_component(c)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")

    result = network.solve(tol=1e-9, max_iter=100)
    assert result.converged

    # Composition passes through unchanged (Pipe doesn't override
    # outlet_fluid()) -- the same HumidAirFluid instance at both ends.
    assert result.node_fluid["src.out"] is fluid
    assert result.node_fluid["pipe.out"] is fluid

    # Not just bookkeeping: Pipe's own pressure-drop residual must have used
    # the humid-air fluid's actual density, not the network's default
    # IdealGasFluid's -- mirrors
    # test_fluid_propagation.py::test_resolve_node_fluid_downstream_component_actually_reads_the_swap.
    P_in, h_in = result.node_P["pipe.in"], result.node_h["pipe.in"]
    rho_humid_air = fluid.density_ph(P_in, h_in)
    rho_dry_air = network.fluid.density_ph(P_in, h_in)
    assert not math.isclose(rho_humid_air, rho_dry_air, rel_tol=1e-3)

    mdot = result.node_mdot["pipe.in"]
    area = math.pi * pipe.D**2 / 4
    v = mdot / (rho_humid_air * area)
    expected_dp = pipe.f * (pipe.L / pipe.D) * (rho_humid_air * v**2 / 2)
    assert math.isclose(
        result.node_P["pipe.in"] - result.node_P["pipe.out"], expected_dp, rel_tol=1e-6
    )
