import pytest

from thermowave.components.heat_transfer import (
    Conduction,
    Convection,
    Radiation,
    ThermalMass,
    heat_loss_watts,
)
from thermowave.core.constants import STEFAN_BOLTZMANN
from thermowave.core.network import NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _FakeState:
    """Minimal stand-in for core.network.NetworkState."""

    def __init__(self, params):
        self._params = params

    def param(self, name: str) -> float:
        return self._params[name]


def test_convection_q_between_two_thermal_masses():
    hot = ThermalMass(name="hot", thermal_capacitance=100.0, T0=400.0)
    cold = ThermalMass(name="cold", thermal_capacitance=100.0, T0=300.0)
    conv = Convection(name="conv1", a=hot, b=cold, h=25.0, A=2.0)
    state = _FakeState({"hot.T": 400.0, "cold.T": 300.0})
    assert conv.Q(state) == pytest.approx(25.0 * 2.0 * (400.0 - 300.0))


def test_convection_q_to_fixed_ambient_temperature():
    hot = ThermalMass(name="hot", thermal_capacitance=100.0, T0=350.0)
    conv = Convection(name="conv1", a=hot, b=288.15, h=10.0, A=1.5)
    state = _FakeState({"hot.T": 350.0})
    assert conv.Q(state) == pytest.approx(10.0 * 1.5 * (350.0 - 288.15))


def test_convection_negative_q_when_b_is_hotter():
    a = ThermalMass(name="a", thermal_capacitance=100.0, T0=300.0)
    b = ThermalMass(name="b", thermal_capacitance=100.0, T0=400.0)
    conv = Convection(name="conv1", a=a, b=b, h=10.0, A=1.0)
    state = _FakeState({"a.T": 300.0, "b.T": 400.0})
    assert conv.Q(state) == pytest.approx(10.0 * 1.0 * (300.0 - 400.0))
    assert conv.Q(state) < 0.0


def test_convection_temperature_source_from_component_port():
    src_P, src_h = 101325.0, AIR.enthalpy_pt(101325.0, 500.0)
    state = NetworkState(
        fluid=AIR,
        node_P={"comp.out": src_P},
        node_h={"comp.out": src_h},
        node_mdot={"comp.out": 1.0},
        params={"ambient_mass.T": 288.15},
    )

    class _FakeFlowComponent:
        name = "comp"

        def ports(self):
            return {"out": "comp.out"}

    comp = _FakeFlowComponent()
    ambient = ThermalMass(name="ambient_mass", thermal_capacitance=100.0, T0=288.15)
    conv = Convection(name="conv1", a=(comp, "out"), b=ambient, h=5.0, A=1.0)

    T_a = AIR.temperature_ph(src_P, src_h)
    assert conv.Q(state) == pytest.approx(5.0 * 1.0 * (T_a - 288.15))


def test_conduction_q_uses_k_a_over_l_conductance():
    a = ThermalMass(name="a", thermal_capacitance=100.0, T0=500.0)
    b = ThermalMass(name="b", thermal_capacitance=100.0, T0=400.0)
    cond = Conduction(name="cond1", a=a, b=b, k=15.0, A=0.05, L=0.2)
    assert cond.UA == pytest.approx(15.0 * 0.05 / 0.2)
    state = _FakeState({"a.T": 500.0, "b.T": 400.0})
    assert cond.Q(state) == pytest.approx(cond.UA * (500.0 - 400.0))


def test_radiation_q_follows_stefan_boltzmann_t4_law():
    a = ThermalMass(name="a", thermal_capacitance=100.0, T0=1000.0)
    b = ThermalMass(name="b", thermal_capacitance=100.0, T0=500.0)
    rad = Radiation(name="rad1", a=a, b=b, emissivity=0.8, A=0.5)
    state = _FakeState({"a.T": 1000.0, "b.T": 500.0})
    expected = 0.8 * 1.0 * STEFAN_BOLTZMANN * 0.5 * (1000.0**4 - 500.0**4)
    assert rad.Q(state) == pytest.approx(expected)


def test_radiation_rejects_out_of_range_emissivity():
    a = ThermalMass(name="a", thermal_capacitance=100.0, T0=1000.0)
    b = ThermalMass(name="b", thermal_capacitance=100.0, T0=500.0)
    with pytest.raises(ValueError, match="emissivity must be in"):
        Radiation(name="rad1", a=a, b=b, emissivity=0.0, A=1.0)
    with pytest.raises(ValueError, match="emissivity must be in"):
        Radiation(name="rad1", a=a, b=b, emissivity=1.5, A=1.0)


def test_radiation_view_factor_scales_q():
    a = ThermalMass(name="a", thermal_capacitance=100.0, T0=1000.0)
    b = ThermalMass(name="b", thermal_capacitance=100.0, T0=500.0)
    rad_full = Radiation(name="rad_full", a=a, b=b, emissivity=1.0, A=1.0, view_factor=1.0)
    rad_half = Radiation(name="rad_half", a=a, b=b, emissivity=1.0, A=1.0, view_factor=0.5)
    state = _FakeState({"a.T": 1000.0, "b.T": 500.0})
    assert rad_half.Q(state) == pytest.approx(0.5 * rad_full.Q(state))


@pytest.mark.parametrize("cls_kwargs", [
    {"cls": Convection, "kwargs": {"h": 10.0, "A": 1.0}},
    {"cls": Conduction, "kwargs": {"k": 10.0, "A": 1.0, "L": 1.0}},
    {"cls": Radiation, "kwargs": {"emissivity": 0.8, "A": 1.0}},
])
def test_heat_path_ports_residuals_and_category(cls_kwargs):
    a = ThermalMass(name="a", thermal_capacitance=100.0, T0=400.0)
    b = ThermalMass(name="b", thermal_capacitance=100.0, T0=300.0)
    path = cls_kwargs["cls"](name="p1", a=a, b=b, **cls_kwargs["kwargs"])
    assert path.ports() == {}
    assert path.residuals(_FakeState({})) == []
    assert path.report_category() == "heat_transfer"


def test_heat_path_report_metrics_exposes_q_and_endpoint_temperatures():
    a = ThermalMass(name="a", thermal_capacitance=100.0, T0=400.0)
    b = ThermalMass(name="b", thermal_capacitance=100.0, T0=300.0)
    conv = Convection(name="conv1", a=a, b=b, h=10.0, A=1.0)
    state = _FakeState({"a.T": 400.0, "b.T": 300.0})
    metrics = conv.report_metrics(state)
    assert metrics["Q [W]"] == pytest.approx(10.0 * 1.0 * (400.0 - 300.0))
    assert metrics["T_a [K]"] == 400.0
    assert metrics["T_b [K]"] == 300.0


def test_heat_loss_watts_none_returns_zero():
    state = _FakeState({})
    assert heat_loss_watts(None, state) == 0.0


def test_heat_loss_watts_delegates_to_path_q():
    a = ThermalMass(name="a", thermal_capacitance=100.0, T0=400.0)
    b = ThermalMass(name="b", thermal_capacitance=100.0, T0=300.0)
    conv = Convection(name="conv1", a=a, b=b, h=10.0, A=1.0)
    state = _FakeState({"a.T": 400.0, "b.T": 300.0})
    assert heat_loss_watts(conv, state) == conv.Q(state)
