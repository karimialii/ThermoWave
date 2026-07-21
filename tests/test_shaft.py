import math

import pytest

from thermowave.components.compressor import Compressor
from thermowave.components.controller import Controller
from thermowave.components.pipe import Pipe
from thermowave.components.sensor import Sensor
from thermowave.components.shaft import Shaft
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)
GAMMA = 1005.0 / (1005.0 - 287.05)


class _FakeState:
    def __init__(self, params):
        self._params = params

    def param(self, name):
        return self._params[name]


class _FakeShaftComponent:
    def __init__(self, name, N=None, power=None):
        self.name = name
        self._N = N
        self._power = power

    def free_parameters(self):
        return {} if self._N is not None else {"N": 50000.0}

    def report_metrics(self, state):
        metrics = {}
        if self._power is not None:
            metrics["power [W]"] = self._power
        return metrics


def _free_pair():
    comp = Compressor(name="c1", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None)
    turb = Turbine(name="t1", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=None)
    return comp, turb


def test_ports_is_empty():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", components=[comp, turb])
    assert shaft.ports() == {}


def test_report_category_is_shaft():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", components=[comp, turb])
    assert shaft.report_category() == "shaft"


def test_raises_if_fewer_than_two_components():
    comp, _turb = _free_pair()
    with pytest.raises(ValueError, match="at least 2"):
        Shaft(name="shaft", components=[comp])


def test_raises_if_any_component_free_param_not_declared_free():
    comp = Compressor(name="c1", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=50000.0)
    _comp2, turb = _free_pair()
    with pytest.raises(ValueError, match="needs at least .* declaring"):
        Shaft(name="shaft", components=[comp, turb])


def test_raises_if_gear_ratios_length_mismatched():
    comp, turb = _free_pair()
    with pytest.raises(ValueError, match="gear_ratios"):
        Shaft(name="shaft", components=[comp, turb], gear_ratios=[1.0, 2.0])


def test_raises_if_signs_length_mismatched():
    comp, turb = _free_pair()
    with pytest.raises(ValueError, match="signs"):
        Shaft(name="shaft", components=[comp, turb], signs=[1.0])


def test_residuals_zero_when_speeds_match_direct_coupling():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", components=[comp, turb])
    state = _FakeState({"c1.N": 50000.0, "t1.N": 50000.0})
    assert math.isclose(shaft.residuals(state)[0], 0.0, abs_tol=1e-9)


def test_residuals_reflect_gear_ratio():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", components=[comp, turb], gear_ratios=[2.0])
    state = _FakeState({"c1.N": 100000.0, "t1.N": 200000.0})
    assert math.isclose(shaft.residuals(state)[0], 0.0, abs_tol=1e-9)


def test_residuals_one_per_follower_for_three_components():
    a = _FakeShaftComponent("a", N=None)
    b = _FakeShaftComponent("b", N=None)
    c = _FakeShaftComponent("c", N=None)
    shaft = Shaft(name="shaft", components=[a, b, c])
    state = _FakeState({"a.N": 1000.0, "b.N": 900.0, "c.N": 1100.0})
    residuals = shaft.residuals(state)
    assert len(residuals) == 2
    assert math.isclose(residuals[0], -100.0)
    assert math.isclose(residuals[1], 100.0)


def test_report_metrics_net_power_uses_signs_and_efficiency():
    a = _FakeShaftComponent("a", N=None, power=1000.0)  # consumer
    b = _FakeShaftComponent("b", N=None, power=1500.0)  # producer
    shaft = Shaft(name="shaft", components=[a, b], signs=[-1.0, 1.0], efficiency=0.9)
    state = _FakeState({"a.N": 5000.0, "b.N": 5000.0})
    metrics = shaft.report_metrics(state)
    assert math.isclose(metrics["power [W]"], 0.9 * (1500.0 - 1000.0))
    assert math.isclose(metrics["eta [-]"], 0.9)
    assert math.isclose(metrics["N [rev/min]"], 5000.0)


def test_report_metrics_includes_inertia():
    a = _FakeShaftComponent("a", N=None)
    b = _FakeShaftComponent("b", N=None)
    shaft = Shaft(name="shaft", components=[a, b], inertia=0.75)
    state = _FakeState({"a.N": 50000.0, "b.N": 50000.0})
    assert math.isclose(shaft.report_metrics(state)["inertia [kg*m^2]"], 0.75)


def test_shaft_end_to_end_locks_compressor_and_turbine_speed():
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    comp = Compressor(name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None)
    heater = Pipe(name="heater", L=1.0, D=0.1, f=0.0, n_elem=1, heat_loss=-100000.0)
    turb = Turbine(name="turb", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=None)
    shaft = Shaft(name="shaft", components=[comp, turb], signs=[-1.0, 1.0], efficiency=0.98)
    snk = Sink(name="snk")

    turb_outlet_sensor = Sensor(name="turb_outlet_sensor")
    target_T = 500.0
    ctrl = Controller(
        name="ctrl",
        sensor=turb_outlet_sensor,
        quantity="T [K]",
        component=comp,
        free_param="N",
        value=target_T,
    )

    network = Network(fluid=AIR)
    for component in (src, comp, heater, turb, shaft, turb_outlet_sensor, ctrl, snk):
        network.add_component(component)

    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", heater, "in")
    network.connect(heater, "out", turb, "in")
    network.connect(turb, "out", turb_outlet_sensor, "tap")
    network.connect(turb, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=400, damping=0.3)
    assert math.isclose(result.params["comp.N"], result.params["turb.N"], rel_tol=1e-6)
    T_out = AIR.temperature_ph(result.node_P["turb.out"], result.node_h["turb.out"])
    assert math.isclose(T_out, target_T, abs_tol=1e-3)


def test_dynamic_shaft_declares_no_differential_state_by_default():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", components=[comp, turb])
    assert shaft.differential_parameters() == {}
    assert shaft.state_derivative(state=None) == {}


def test_dynamic_shaft_raises_if_inertia_not_positive():
    comp, turb = _free_pair()
    with pytest.raises(ValueError, match="inertia"):
        Shaft(name="shaft", components=[comp, turb], dynamic=True, inertia=0.0)


def test_dynamic_shaft_gear_ratios_need_one_entry_per_component_not_follower():
    comp, turb = _free_pair()
    with pytest.raises(ValueError, match="gear_ratios"):
        Shaft(
            name="shaft", components=[comp, turb], dynamic=True, inertia=0.05,
            gear_ratios=[1.0],  # would be valid for dynamic=False, not dynamic=True
        )


def test_dynamic_shaft_declares_differential_parameter_seeded_by_n0():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", components=[comp, turb], dynamic=True, inertia=0.05, N0=42000.0)
    assert shaft.differential_parameters() == {"N": 42000.0}


def test_dynamic_shaft_residuals_couple_every_component_to_shaft_speed():
    a = _FakeShaftComponent("a", N=None)
    b = _FakeShaftComponent("b", N=None)
    c = _FakeShaftComponent("c", N=None)
    shaft = Shaft(
        name="shaft", components=[a, b, c], dynamic=True, inertia=0.05,
        gear_ratios=[1.0, 2.0, 0.5],
    )
    state = _FakeState({"a.N": 1000.0, "b.N": 2000.0, "c.N": 500.0, "shaft.N": 1000.0})
    residuals = shaft.residuals(state)
    assert len(residuals) == 3
    for residual in residuals:
        assert math.isclose(residual, 0.0, abs_tol=1e-9)


def test_dynamic_shaft_state_derivative_is_zero_when_powers_balance():
    a = _FakeShaftComponent("a", N=None, power=1000.0)  # consumer
    b = _FakeShaftComponent("b", N=None, power=1000.0)  # producer, exactly balances
    shaft = Shaft(
        name="shaft", components=[a, b], signs=[-1.0, 1.0],
        dynamic=True, inertia=0.05,
    )
    state = _FakeState({"a.N": 5000.0, "b.N": 5000.0, "shaft.N": 5000.0})
    assert math.isclose(shaft.state_derivative(state)["N"], 0.0, abs_tol=1e-9)


def test_dynamic_shaft_state_derivative_positive_when_net_power_positive():
    a = _FakeShaftComponent("a", N=None, power=1000.0)  # consumer
    b = _FakeShaftComponent("b", N=None, power=1500.0)  # producer, net positive
    shaft = Shaft(
        name="shaft", components=[a, b], signs=[-1.0, 1.0],
        dynamic=True, inertia=0.05,
    )
    state = _FakeState({"a.N": 5000.0, "b.N": 5000.0, "shaft.N": 5000.0})
    assert shaft.state_derivative(state)["N"] > 0.0


def test_dynamic_shaft_report_metrics_reads_speed_from_differential_state():
    a = _FakeShaftComponent("a", N=None)
    b = _FakeShaftComponent("b", N=None)
    shaft = Shaft(name="shaft", components=[a, b], dynamic=True, inertia=0.05)
    state = _FakeState({"a.N": 5000.0, "b.N": 5000.0, "shaft.N": 4321.0})
    assert math.isclose(shaft.report_metrics(state)["N [rev/min]"], 4321.0)
