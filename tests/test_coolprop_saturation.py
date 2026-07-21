import math

import pytest

pytest.importorskip("CoolProp")

from thermowave.fluids.real_fluid import CoolPropFluid  # noqa: E402
from thermowave.fluids.two_phase import (  # noqa: E402
    require_entropy,
    require_two_phase,
    supports_entropy,
    supports_two_phase,
)
from thermowave.fluids.ideal_gas import IdealGasFluid  # noqa: E402

WATER = CoolPropFluid(name="Water")
P_ATM = 101325.0


def test_saturation_temperature_of_water_at_one_atm():
    assert math.isclose(WATER.saturation_temperature(P_ATM), 373.124, abs_tol=0.01)


def test_saturation_pressure_round_trips_saturation_temperature():
    T_sat = WATER.saturation_temperature(P_ATM)
    assert math.isclose(WATER.saturation_pressure(T_sat), P_ATM, rel_tol=1e-4)


def test_saturated_vapor_enthalpy_exceeds_saturated_liquid():
    h_f = WATER.saturated_liquid_enthalpy(P_ATM)
    h_g = WATER.saturated_vapor_enthalpy(P_ATM)
    assert h_g > h_f
    # Latent heat of vaporization of water at 1 atm ~ 2.26 MJ/kg.
    assert math.isclose(h_g - h_f, 2.256e6, rel_tol=1e-2)


def test_enthalpy_pq_midpoint_is_mean_of_hf_and_hg():
    h_f = WATER.saturated_liquid_enthalpy(P_ATM)
    h_g = WATER.saturated_vapor_enthalpy(P_ATM)
    assert math.isclose(WATER.enthalpy_pq(P_ATM, 0.5), 0.5 * (h_f + h_g), rel_tol=1e-9)


def test_quality_ph_round_trips_enthalpy_pq():
    for x in (0.0, 0.25, 0.5, 0.9, 1.0):
        h = WATER.enthalpy_pq(P_ATM, x)
        assert math.isclose(WATER.quality_ph(P_ATM, h), x, abs_tol=1e-6)


def test_quality_ph_is_minus_one_outside_the_dome():
    # CoolProp returns -1.0 for BOTH subcooled and superheated single-phase.
    h_subcooled = WATER.enthalpy_pt(P_ATM, 300.0)
    h_superheated = WATER.enthalpy_pt(P_ATM, 500.0)
    assert WATER.quality_ph(P_ATM, h_subcooled) == -1.0
    assert WATER.quality_ph(P_ATM, h_superheated) == -1.0


def test_enthalpy_ps_round_trips_entropy_ph():
    h = 2.0e6  # J/kg, inside the two-phase region at 1 atm
    s = WATER.entropy_ph(P_ATM, h)
    assert math.isclose(WATER.enthalpy_ps(P_ATM, s), h, rel_tol=1e-9)


def test_supports_two_phase_true_for_coolprop_false_for_ideal_gas():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    assert supports_two_phase(WATER)
    assert supports_entropy(WATER)
    assert not supports_two_phase(air)
    assert not supports_entropy(air)


def test_require_two_phase_raises_for_ideal_gas():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    with pytest.raises(ValueError, match="two-phase-capable"):
        require_two_phase(air, "SimpleEvaporator")
    require_two_phase(WATER, "SimpleEvaporator")  # no raise


def test_require_entropy_raises_for_ideal_gas():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    with pytest.raises(ValueError, match="entropy"):
        require_entropy(air, "Pump")
    require_entropy(WATER, "Pump")  # no raise
