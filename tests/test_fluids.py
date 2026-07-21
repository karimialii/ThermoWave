import math

import pytest

from thermowave.fluids.ideal_gas import IdealGasFluid

AIR_R = 287.05  # J/(kg*K)
AIR_CP = 1005.0  # J/(kg*K)


def test_ideal_gas_enthalpy_pt_matches_hand_calculation():
    air = IdealGasFluid(name="air", R=AIR_R, cp=AIR_CP)
    h = air.enthalpy_pt(P=101325.0, T=300.0)
    assert math.isclose(h, AIR_CP * 300.0, rel_tol=1e-9)


def test_ideal_gas_temperature_ph_matches_hand_calculation():
    air = IdealGasFluid(name="air", R=AIR_R, cp=AIR_CP)
    h = AIR_CP * 300.0
    T = air.temperature_ph(P=101325.0, h=h)
    assert math.isclose(T, 300.0, rel_tol=1e-9)


def test_ideal_gas_density_ph_matches_pv_equals_mrt():
    air = IdealGasFluid(name="air", R=AIR_R, cp=AIR_CP)
    P = 101325.0
    T = 300.0
    h = AIR_CP * T
    rho = air.density_ph(P=P, h=h)
    expected_rho = P / (AIR_R * T)
    assert math.isclose(rho, expected_rho, rel_tol=1e-9)


def test_ideal_gas_cp_is_constant():
    air = IdealGasFluid(name="air", R=AIR_R, cp=AIR_CP)
    assert air.cp(P=101325.0, T=300.0) == AIR_CP
    assert air.cp(P=200000.0, T=500.0) == AIR_CP


def test_ideal_gas_cv_matches_cp_minus_r():
    air = IdealGasFluid(name="air", R=AIR_R, cp=AIR_CP)
    assert math.isclose(air.cv(P=101325.0, T=300.0), AIR_CP - AIR_R, rel_tol=1e-12)


def test_ideal_gas_gamma_matches_cp_over_cv():
    air = IdealGasFluid(name="air", R=AIR_R, cp=AIR_CP)
    expected_gamma = AIR_CP / (AIR_CP - AIR_R)
    assert math.isclose(air.gamma(P=101325.0, T=300.0), expected_gamma, rel_tol=1e-12)
    assert math.isclose(expected_gamma, 1.4, rel_tol=1e-2)  # ~1.4 for air, sanity check


def test_ideal_gas_name_attribute():
    air = IdealGasFluid(name="air", R=AIR_R, cp=AIR_CP)
    assert air.name == "air"


from thermowave.fluids.ideal_gas_mixture import IdealGasMixtureFluid

_R_UNIVERSAL = 8.314462618


def test_ideal_gas_mixture_matches_plain_ideal_gas_for_single_species():
    n2 = IdealGasMixtureFluid(name="n2", composition={"N2": 1.0})
    M_n2, cp_molar_n2 = IdealGasMixtureFluid.SPECIES["N2"]
    expected_R = _R_UNIVERSAL / M_n2
    expected_cp = cp_molar_n2 / M_n2
    assert math.isclose(n2.R, expected_R, rel_tol=1e-9)
    assert math.isclose(n2._cp, expected_cp, rel_tol=1e-9)
    h = n2.enthalpy_pt(P=101325.0, T=400.0)
    assert math.isclose(h, expected_cp * 400.0, rel_tol=1e-9)
    T_back = n2.temperature_ph(P=101325.0, h=h)
    assert math.isclose(T_back, 400.0, rel_tol=1e-9)


def test_ideal_gas_mixture_r_and_cp_are_mass_fraction_weighted():
    air = IdealGasMixtureFluid(name="air", composition={"O2": 0.233, "N2": 0.767})
    M_o2, cp_o2 = IdealGasMixtureFluid.SPECIES["O2"]
    M_n2, cp_n2 = IdealGasMixtureFluid.SPECIES["N2"]
    expected_R = 0.233 * (_R_UNIVERSAL / M_o2) + 0.767 * (_R_UNIVERSAL / M_n2)
    expected_cp = 0.233 * (cp_o2 / M_o2) + 0.767 * (cp_n2 / M_n2)
    assert math.isclose(air.R, expected_R, rel_tol=1e-9)
    assert math.isclose(air._cp, expected_cp, rel_tol=1e-9)


def test_ideal_gas_mixture_density_ph_matches_pv_equals_mrt():
    air = IdealGasMixtureFluid(name="air", composition={"O2": 0.233, "N2": 0.767})
    P, T = 200000.0, 350.0
    h = air.enthalpy_pt(P, T)
    rho = air.density_ph(P, h)
    assert math.isclose(rho, P / (air.R * T), rel_tol=1e-9)


def test_ideal_gas_mixture_raises_if_composition_does_not_sum_to_one():
    with pytest.raises(ValueError, match="sum to"):
        IdealGasMixtureFluid(name="bad", composition={"O2": 0.5, "N2": 0.3})


def test_ideal_gas_mixture_raises_on_unknown_species():
    with pytest.raises(ValueError, match="unknown species"):
        IdealGasMixtureFluid(name="bad", composition={"Xenon": 1.0})


def test_ideal_gas_mixture_extra_species_can_add_a_custom_entry():
    fluid = IdealGasMixtureFluid(
        name="propane_air", composition={"C3H8": 0.05, "N2": 0.95},
        extra_species={"C3H8": (0.044097, 73.6)},
    )
    M, cp_molar = 0.044097, 73.6
    M_n2, cp_n2 = IdealGasMixtureFluid.SPECIES["N2"]
    expected_R = 0.05 * (_R_UNIVERSAL / M) + 0.95 * (_R_UNIVERSAL / M_n2)
    assert math.isclose(fluid.R, expected_R, rel_tol=1e-9)


pytest.importorskip("cantera")

from thermowave.fluids.cantera_fluid import CanteraFluid  # noqa: E402


def test_cantera_fluid_enthalpy_pt_matches_direct_cantera_call():
    import cantera as ct

    air = CanteraFluid(name="air", composition="O2:0.21, N2:0.79")
    P, T = 150000.0, 450.0
    h = air.enthalpy_pt(P, T)

    gas = ct.Solution("gri30.yaml")
    gas.TPX = T, P, "O2:0.21, N2:0.79"
    assert math.isclose(h, gas.enthalpy_mass, rel_tol=1e-9)


def test_cantera_fluid_round_trip_temperature_enthalpy():
    air = CanteraFluid(name="air", composition="O2:0.21, N2:0.79")
    P, T = 150000.0, 450.0
    h = air.enthalpy_pt(P, T)
    T_back = air.temperature_ph(P, h)
    assert math.isclose(T_back, T, rel_tol=1e-6)


def test_cantera_fluid_cp_is_positive_and_reasonable_for_air():
    air = CanteraFluid(name="air", composition="O2:0.21, N2:0.79")
    cp = air.cp(P=101325.0, T=300.0)
    assert 900.0 < cp < 1100.0  # air's cp is ~1005 J/(kg*K) near 300 K


def test_cantera_fluid_cv_matches_direct_cantera_call():
    import cantera as ct

    air = CanteraFluid(name="air", composition="O2:0.21, N2:0.79")
    P, T = 150000.0, 450.0
    cv = air.cv(P, T)

    gas = ct.Solution("gri30.yaml")
    gas.TPX = T, P, "O2:0.21, N2:0.79"
    assert math.isclose(cv, gas.cv_mass, rel_tol=1e-9)


def test_cantera_fluid_gamma_is_reasonable_for_air():
    air = CanteraFluid(name="air", composition="O2:0.21, N2:0.79")
    gamma = air.gamma(P=101325.0, T=300.0)
    assert 1.35 < gamma < 1.45  # air's gamma is ~1.4 near 300 K


def test_cantera_fluid_density_ph_matches_ideal_gas_at_low_pressure():
    air = CanteraFluid(name="air", composition="O2:0.21, N2:0.79")
    P, T = 101325.0, 300.0
    h = air.enthalpy_pt(P, T)
    rho = air.density_ph(P, h)
    R_air = 287.05
    assert math.isclose(rho, P / (R_air * T), rel_tol=1e-2)


def test_cantera_fluid_mass_basis_accepts_mass_fraction_composition():
    air = CanteraFluid(name="air", composition="O2:0.232, N2:0.768", basis="mass")
    h = air.enthalpy_pt(P=101325.0, T=300.0)
    assert math.isfinite(h)


def test_cantera_fluid_rejects_invalid_basis():
    with pytest.raises(ValueError, match="basis"):
        CanteraFluid(name="air", composition="O2:0.21, N2:0.79", basis="volume")


def test_cantera_fluid_import_error_message_when_cantera_missing(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "cantera" or name.startswith("cantera"):
            raise ImportError("simulated missing cantera")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="cantera"):
        CanteraFluid(name="air", composition="O2:0.21, N2:0.79")


pytest.importorskip("CoolProp")

from thermowave.core.exceptions import FluidRangeError  # noqa: E402
from thermowave.fluids.real_fluid import CoolPropFluid  # noqa: E402


def test_coolprop_fluid_enthalpy_pt_matches_direct_propsi_call():
    from CoolProp.CoolProp import PropsSI

    water = CoolPropFluid(name="Water")
    P = 1.0e5  # 1 bar
    T = 473.15  # 200 C, superheated steam at 1 bar
    h = water.enthalpy_pt(P=P, T=T)
    expected = PropsSI("H", "P", P, "T", T, "Water")
    assert h == pytest.approx(expected, rel=1e-9)


def test_coolprop_fluid_round_trip_temperature_enthalpy():
    water = CoolPropFluid(name="Water")
    P = 1.0e5
    T = 473.15
    h = water.enthalpy_pt(P=P, T=T)
    T_back = water.temperature_ph(P=P, h=h)
    assert T_back == pytest.approx(T, rel=1e-6)


def test_coolprop_fluid_cv_matches_direct_propsi_call():
    from CoolProp.CoolProp import PropsSI

    water = CoolPropFluid(name="Water")
    P = 1.0e5
    T = 473.15  # superheated steam at 1 bar, away from the two-phase dome
    cv = water.cv(P=P, T=T)
    expected = PropsSI("O", "P", P, "T", T, "Water")
    assert cv == pytest.approx(expected, rel=1e-9)


def test_coolprop_fluid_density_ph_matches_direct_propsi_call():
    from CoolProp.CoolProp import PropsSI

    water = CoolPropFluid(name="Water")
    P = 1.0e5
    T = 473.15
    h = water.enthalpy_pt(P=P, T=T)
    rho = water.density_ph(P=P, h=h)
    expected = PropsSI("D", "P", P, "H", h, "Water")
    assert rho == pytest.approx(expected, rel=1e-9)


def test_coolprop_fluid_clamps_pressure_below_minimum():
    water = CoolPropFluid(name="Water", P_min=1.0e5, P_max=1.0e7)
    assert water._clamp_pressure(1.0) == 1.0e5
    assert water._clamp_pressure(1.0e9) == 1.0e7
    assert water._clamp_pressure(5.0e5) == 5.0e5


def test_coolprop_fluid_raises_fluid_range_error_on_invalid_property_call():
    water = CoolPropFluid(name="Water")
    with pytest.raises(FluidRangeError):
        water.density_ph(P=1.0e5, h=-1.0e9)


def test_coolprop_fluid_import_error_message_when_coolprop_missing(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "CoolProp.CoolProp" or name.startswith("CoolProp"):
            raise ImportError("simulated missing CoolProp")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="coolprop"):
        CoolPropFluid(name="Water")
