import math

import pytest

from thermowave.components.compressor import Compressor
from thermowave.components.controller import Controller
from thermowave.components.sensor import Sensor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _FakeState:
    def __init__(self, fluid, mdot, node_values):
        self.fluid = fluid
        self._mdot = mdot
        self._node_values = node_values

    def node(self, name):
        return self._node_values[name]

    def fluid_at(self, name):
        return self.fluid

    def mdot(self, name):
        return self._mdot[name]


def test_sensor_ports_returns_single_tap_port_derived_from_name():
    sensor = Sensor(name="s1")
    assert sensor.ports() == {"tap": "s1.tap"}


def test_sensor_residuals_is_empty():
    sensor = Sensor(name="s1")
    assert sensor.residuals(state=None) == []


def test_sensor_report_metrics_reads_node_p_t_h_and_mdot():
    air = AIR
    P, T = 300000.0, 400.0
    h = air.enthalpy_pt(P, T)
    sensor = Sensor(name="s1")
    state = _FakeState(fluid=air, mdot={"s1.tap": 0.5}, node_values={"s1.tap": (P, h)})
    metrics = sensor.report_metrics(state)
    assert math.isclose(metrics["P [Pa]"], P)
    assert math.isclose(metrics["T [K]"], T)
    assert math.isclose(metrics["h [J/kg]"], h)
    assert math.isclose(metrics["mdot [kg/s]"], 0.5)


def test_sensor_report_metrics_omits_mdot_when_not_available():
    air = AIR
    P, T = 300000.0, 400.0
    h = air.enthalpy_pt(P, T)
    sensor = Sensor(name="s1")

    class _NoMdotState(_FakeState):
        def mdot(self, name):
            raise KeyError(name)

    state = _NoMdotState(fluid=air, mdot={}, node_values={"s1.tap": (P, h)})
    metrics = sensor.report_metrics(state)
    assert "mdot [kg/s]" not in metrics


def test_controller_raises_if_free_param_not_declared_free():
    comp = Compressor(name="c1", map_path="tests/fixtures/simple_compressor_map.cop", gamma=1.4, N=50000.0)
    sensor = Sensor(name="s1")
    with pytest.raises(ValueError, match="doesn't currently declare"):
        Controller(
            name="ctrl1",
            sensor=sensor,
            quantity="T [K]",
            component=comp,
            free_param="N",
            value=400.0,
        )


def test_controller_residual_is_measured_minus_target():
    comp = Compressor(name="c1", map_path="tests/fixtures/simple_compressor_map.cop", gamma=1.4, N=None)
    sensor = Sensor(name="s1")
    ctrl = Controller(
        name="ctrl1", sensor=sensor, quantity="T [K]", component=comp, free_param="N", value=400.0
    )
    air = AIR
    P, T = 300000.0, 450.0
    h = air.enthalpy_pt(P, T)
    state = _FakeState(fluid=air, mdot={}, node_values={"s1.tap": (P, h)})
    assert math.isclose(ctrl.residuals(state)[0], 50.0, abs_tol=1e-6)


def test_controller_ports_is_empty():
    comp = Compressor(name="c1", map_path="tests/fixtures/simple_compressor_map.cop", gamma=1.4, N=None)
    sensor = Sensor(name="s1")
    ctrl = Controller(
        name="ctrl1", sensor=sensor, quantity="T [K]", component=comp, free_param="N", value=400.0
    )
    assert ctrl.ports() == {}


def test_controller_drives_compressor_outlet_temperature_end_to_end():
    gamma = 1005.0 / (1005.0 - 287.05)
    air = AIR

    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    comp = Compressor(name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=gamma, N=None)
    snk = Sink(name="snk")
    sensor = Sensor(name="outlet_sensor")

    target_T = 420.0
    ctrl = Controller(
        name="ctrl",
        sensor=sensor,
        quantity="T [K]",
        component=comp,
        free_param="N",
        value=target_T,
    )

    network = Network(fluid=air)
    for component in (src, comp, sensor, ctrl, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", sensor, "tap")
    network.connect(comp, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=200, damping=0.5)
    P_out, h_out = result.node_P["comp.out"], result.node_h["comp.out"]
    T_out = air.temperature_ph(P_out, h_out)
    assert math.isclose(T_out, target_T, abs_tol=1e-4)
