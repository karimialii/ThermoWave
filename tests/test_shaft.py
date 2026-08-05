import math
import warnings

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
    def __init__(self, node_N=None, params=None):
        self._node_N = node_N or {}
        self._params = params or {}

    def N(self, name):
        return self._node_N[name]

    def param(self, name):
        return self._params[name]

    def signal(self, name):
        return self._node_N[name]  # unused key space reused for simplicity


class _FakeShaftComponent:
    def __init__(self, name, N=None, power=None, has_shaft=True):
        self.name = name
        self._N = N
        self._power = power
        self._has_shaft = has_shaft
        self._shaft_port = f"{name}.shaft"
        self._power_port = f"{name}.power"

    def ports(self):
        return {}

    def mechanical_ports(self):
        return {"shaft": self._shaft_port} if self._has_shaft else {}

    def signal_ports(self):
        return {"power": self._power_port}

    def free_mechanical_ports(self):
        if not self._has_shaft or self._N is not None:
            return {}
        return {self._shaft_port: 50000.0}

    def provided_signal_values(self, state):
        return {} if self._power is None else {self._power_port: self._power}

    def shaft_sign(self):
        return 1.0


def _free_pair():
    comp = Compressor(name="c1", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None)
    turb = Turbine(name="t1", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=None)
    return comp, turb


def _state_for(shaft, N_by_member_index, signal_by_member_index=None):
    """Build a _FakeState keyed by this shaft's own local mechanical/signal
    port ids, in member order -- mirrors how, after a real connect(), the
    shaft's own port id and the member's port id resolve to the same value.
    """
    node_N = {
        list(shaft._mech_ports.values())[i]: N for i, N in enumerate(N_by_member_index)
    }
    if signal_by_member_index is not None:
        for i, power in enumerate(signal_by_member_index):
            node_N[list(shaft._signal_port_names.values())[i]] = power
    return _FakeState(node_N=node_N)


def test_ports_is_empty():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", members=[comp, turb])
    assert shaft.ports() == {}


def test_report_category_is_shaft():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", members=[comp, turb])
    assert shaft.report_category() == "shaft"


def test_raises_if_fewer_than_two_components():
    comp, _turb = _free_pair()
    with pytest.raises(ValueError, match="at least 2"):
        Shaft(name="shaft", members=[comp])


def test_raises_if_any_component_free_param_not_declared_free():
    comp = _FakeShaftComponent("a", has_shaft=True)
    other = _FakeShaftComponent("b", has_shaft=False)
    with pytest.raises(ValueError, match="needs at least .* member"):
        Shaft(name="shaft", members=[comp, other])


def test_signs_are_inferred_from_member_shaft_sign_when_omitted():
    comp, turb = _free_pair()  # real Compressor/Turbine: -1.0 / +1.0
    shaft = Shaft(name="shaft", members=[turb, comp])
    assert shaft.signs == [1.0, -1.0]


def test_explicit_signs_still_override_the_inferred_default():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", members=[turb, comp], signs=[-1.0, -1.0])
    assert shaft.signs == [-1.0, -1.0]


def test_signs_omitted_raises_a_clear_error_for_a_member_with_no_default():
    # _FakeShaftComponent's shaft_sign() defaults to 1.0, so build one that
    # explicitly has none (mirrors a real custom component type that never
    # overrode shaft_sign()) to prove the error names it, rather than
    # crashing or silently guessing.
    class _NoDefaultSign(_FakeShaftComponent):
        def shaft_sign(self):
            return None

    comp, _turb = _free_pair()
    no_default = _NoDefaultSign("mystery")
    with pytest.raises(ValueError, match=r"can't infer signs for \['mystery'\]"):
        Shaft(name="shaft", members=[comp, no_default])


def test_raises_if_gear_ratios_length_mismatched():
    comp, turb = _free_pair()
    with pytest.raises(ValueError, match="gear_ratios"):
        Shaft(name="shaft", members=[comp, turb], gear_ratios=[1.0, 2.0])


def test_raises_if_signs_length_mismatched():
    comp, turb = _free_pair()
    with pytest.raises(ValueError, match="signs"):
        Shaft(name="shaft", members=[comp, turb], signs=[1.0])


def test_residuals_zero_when_speeds_match_direct_coupling():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", members=[comp, turb])
    state = _state_for(shaft, [50000.0, 50000.0])
    assert math.isclose(shaft.residuals(state)[0], 0.0, abs_tol=1e-9)


def test_residuals_reflect_gear_ratio():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", members=[comp, turb], gear_ratios=[2.0])
    state = _state_for(shaft, [100000.0, 200000.0])
    assert math.isclose(shaft.residuals(state)[0], 0.0, abs_tol=1e-9)


def test_residuals_one_per_follower_for_three_components():
    a = _FakeShaftComponent("a", N=None)
    b = _FakeShaftComponent("b", N=None)
    c = _FakeShaftComponent("c", N=None)
    shaft = Shaft(name="shaft", members=[a, b, c])
    state = _state_for(shaft, [1000.0, 900.0, 1100.0])
    residuals = shaft.residuals(state)
    assert len(residuals) == 2
    assert math.isclose(residuals[0], -100.0)
    assert math.isclose(residuals[1], 100.0)


def test_report_metrics_net_power_uses_signs_and_efficiency():
    a = _FakeShaftComponent("a", N=None, power=1000.0)  # consumer
    b = _FakeShaftComponent("b", N=None, power=1500.0)  # producer
    shaft = Shaft(name="shaft", members=[a, b], signs=[-1.0, 1.0], efficiency=0.9)
    state = _state_for(shaft, [5000.0, 5000.0], signal_by_member_index=[1000.0, 1500.0])
    metrics = shaft.report_metrics(state)
    assert math.isclose(metrics["power [W]"], 0.9 * (1500.0 - 1000.0))
    assert math.isclose(metrics["eta [-]"], 0.9)
    assert math.isclose(metrics["N [rev/min]"], 5000.0)


def test_report_metrics_includes_inertia():
    a = _FakeShaftComponent("a", N=None)
    b = _FakeShaftComponent("b", N=None)
    shaft = Shaft(name="shaft", members=[a, b], inertia=0.75)
    state = _state_for(shaft, [50000.0, 50000.0], signal_by_member_index=[0.0, 0.0])
    assert math.isclose(shaft.report_metrics(state)["inertia [kg*m^2]"], 0.75)


def test_shaft_end_to_end_locks_compressor_and_turbine_speed():
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    comp = Compressor(name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None)
    heater = Pipe(name="heater", L=1.0, D=0.1, f=0.0, n_elem=1, heat_loss=-100000.0)
    turb = Turbine(name="turb", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=None)
    shaft = Shaft(name="shaft", members=[comp, turb], signs=[-1.0, 1.0], efficiency=0.98)
    snk = Sink(name="snk")

    turb_outlet_sensor = Sensor(name="turb_outlet_sensor")
    target_T = 500.0
    ctrl = Controller(
        name="ctrl",
        sensor=turb_outlet_sensor,
        quantity="T [K]",
        component=comp,
        free_param="shaft",
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
    # Connecting each member's mechanical "shaft" port pulls its "power"
    # signal port in automatically (Shaft.on_connected()) -- no separate
    # kind="signal" connect() needed for comp/turb.
    network.connect(shaft, "m0", comp, "shaft", kind="mechanical")
    network.connect(shaft, "m1", turb, "shaft", kind="mechanical")

    result = network.solve(tol=1e-8, max_iter=400, damping=0.3)
    assert math.isclose(result.node_N["comp.shaft"], result.node_N["turb.shaft"], rel_tol=1e-6)
    T_out = AIR.temperature_ph(result.node_P["turb.out"], result.node_h["turb.out"])
    assert math.isclose(T_out, target_T, abs_tol=1e-3)


def test_mechanical_connect_alone_wires_power_automatically():
    # The whole point of on_connected(): one connect(kind="mechanical")
    # call is enough, no separate signal connect(), and the network still
    # solves and reports power correctly.
    comp = Compressor(name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None)
    heater = Pipe(name="heater", L=1.0, D=0.1, f=0.0, n_elem=1, heat_loss=-100000.0)
    turb = Turbine(name="turb", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=None)
    shaft = Shaft(name="shaft", members=[comp, turb], signs=[-1.0, 1.0])
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    snk = Sink(name="snk")

    turb_outlet_sensor = Sensor(name="turb_outlet_sensor")
    ctrl = Controller(
        name="ctrl", sensor=turb_outlet_sensor, quantity="T [K]", component=comp,
        free_param="shaft", value=500.0,
    )

    network = Network(fluid=AIR)
    for component in (src, comp, heater, turb, shaft, turb_outlet_sensor, ctrl, snk):
        network.add_component(component)

    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", heater, "in")
    network.connect(heater, "out", turb, "in")
    network.connect(turb, "out", turb_outlet_sensor, "tap")
    network.connect(turb, "out", snk, "in")
    network.connect(shaft, "m0", comp, "shaft", kind="mechanical")
    network.connect(shaft, "m1", turb, "shaft", kind="mechanical")

    kinds = {(c.from_component.name, c.from_port, c.kind) for c in network.connections}
    assert ("shaft", "p0", "signal") in kinds
    assert ("shaft", "p1", "signal") in kinds

    result = network.solve(tol=1e-8, max_iter=400, damping=0.3)
    metrics = shaft.report_metrics(result.state())
    assert math.isfinite(metrics["power [W]"])


def test_autowire_makes_the_same_connections_as_the_manual_pattern():
    # Same network as test_shaft_end_to_end_locks_compressor_and_turbine_speed
    # above, but using shaft.autowire(network) in place of the four explicit
    # network.connect(shaft, "m0"/"m1"/"p0"/"p1", ..., kind=...) calls --
    # should converge to the identical result.
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    comp = Compressor(name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None)
    heater = Pipe(name="heater", L=1.0, D=0.1, f=0.0, n_elem=1, heat_loss=-100000.0)
    turb = Turbine(name="turb", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=None)
    shaft = Shaft(name="shaft", members=[comp, turb], signs=[-1.0, 1.0], efficiency=0.98)
    snk = Sink(name="snk")

    turb_outlet_sensor = Sensor(name="turb_outlet_sensor")
    target_T = 500.0
    ctrl = Controller(
        name="ctrl",
        sensor=turb_outlet_sensor,
        quantity="T [K]",
        component=comp,
        free_param="shaft",
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
    shaft.autowire(network)

    result = network.solve(tol=1e-8, max_iter=400, damping=0.3)
    assert math.isclose(result.node_N["comp.shaft"], result.node_N["turb.shaft"], rel_tol=1e-6)
    T_out = AIR.temperature_ph(result.node_P["turb.out"], result.node_h["turb.out"])
    assert math.isclose(T_out, target_T, abs_tol=1e-3)


def test_autowire_skips_mechanical_port_for_torque_only_members():
    # ShaftLoad has no "shaft" mechanical port -- autowire must not try to
    # connect one for it (that would raise NetworkTopologyError), only the
    # signal port, exactly like the manual pattern in the class docstring.
    from thermowave.components.shaft_load import ShaftLoad

    comp = Compressor(name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=50000.0)
    turb = Turbine(name="turb", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=50000.0)
    load = ShaftLoad(name="load", power=1000.0)
    shaft = Shaft(name="shaft", members=[comp, turb, load], signs=[-1.0, 1.0, -1.0])

    network = Network(fluid=AIR)
    for component in (comp, turb, load, shaft):
        network.add_component(component)
    shaft.autowire(network)

    kinds = {(c.from_component.name, c.from_port, c.kind) for c in network.connections}
    assert ("shaft", "m0", "mechanical") in kinds
    assert ("shaft", "m1", "mechanical") in kinds
    assert ("shaft", "p0", "signal") in kinds
    assert ("shaft", "p1", "signal") in kinds
    assert ("shaft", "p2", "signal") in kinds
    assert not any(port.startswith("m2") for _, port, _ in kinds)


def test_dynamic_shaft_declares_no_differential_state_by_default():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", members=[comp, turb])
    assert shaft.differential_parameters() == {}
    assert shaft.state_derivative(state=None) == {}


def test_dynamic_shaft_raises_if_inertia_not_positive():
    comp, turb = _free_pair()
    with pytest.raises(ValueError, match="inertia"):
        Shaft(name="shaft", members=[comp, turb], dynamic=True, inertia=0.0)


def test_dynamic_shaft_gear_ratios_need_one_entry_per_component_not_follower():
    comp, turb = _free_pair()
    with pytest.raises(ValueError, match="gear_ratios"):
        Shaft(
            name="shaft", members=[comp, turb], dynamic=True, inertia=0.05,
            gear_ratios=[1.0],  # would be valid for dynamic=False, not dynamic=True
        )


def test_dynamic_shaft_declares_differential_parameter_seeded_by_n0():
    comp, turb = _free_pair()
    shaft = Shaft(name="shaft", members=[comp, turb], dynamic=True, inertia=0.05, N0=42000.0)
    assert shaft.differential_parameters() == {"N": 42000.0}


def test_dynamic_shaft_residuals_couple_every_component_to_shaft_speed():
    a = _FakeShaftComponent("a", N=None)
    b = _FakeShaftComponent("b", N=None)
    c = _FakeShaftComponent("c", N=None)
    shaft = Shaft(
        name="shaft", members=[a, b, c], dynamic=True, inertia=0.05,
        gear_ratios=[1.0, 2.0, 0.5],
    )
    state = _state_for(shaft, [1000.0, 2000.0, 500.0])
    state._params["shaft.N"] = 1000.0
    residuals = shaft.residuals(state)
    assert len(residuals) == 3
    for residual in residuals:
        assert math.isclose(residual, 0.0, abs_tol=1e-9)


def test_dynamic_shaft_state_derivative_is_zero_when_powers_balance():
    a = _FakeShaftComponent("a", N=None, power=1000.0)  # consumer
    b = _FakeShaftComponent("b", N=None, power=1000.0)  # producer, exactly balances
    shaft = Shaft(
        name="shaft", members=[a, b], signs=[-1.0, 1.0],
        dynamic=True, inertia=0.05,
    )
    state = _state_for(shaft, [5000.0, 5000.0], signal_by_member_index=[1000.0, 1000.0])
    state._params["shaft.N"] = 5000.0
    assert math.isclose(shaft.state_derivative(state)["N"], 0.0, abs_tol=1e-9)


def test_dynamic_shaft_state_derivative_positive_when_net_power_positive():
    a = _FakeShaftComponent("a", N=None, power=1000.0)  # consumer
    b = _FakeShaftComponent("b", N=None, power=1500.0)  # producer, net positive
    shaft = Shaft(
        name="shaft", members=[a, b], signs=[-1.0, 1.0],
        dynamic=True, inertia=0.05,
    )
    state = _state_for(shaft, [5000.0, 5000.0], signal_by_member_index=[1000.0, 1500.0])
    state._params["shaft.N"] = 5000.0
    assert shaft.state_derivative(state)["N"] > 0.0


def test_dynamic_shaft_report_metrics_reads_speed_from_differential_state():
    a = _FakeShaftComponent("a", N=None)
    b = _FakeShaftComponent("b", N=None)
    shaft = Shaft(name="shaft", members=[a, b], dynamic=True, inertia=0.05)
    state = _state_for(shaft, [5000.0, 5000.0], signal_by_member_index=[0.0, 0.0])
    state._params["shaft.N"] = 4321.0
    assert math.isclose(shaft.report_metrics(state)["N [rev/min]"], 4321.0)


def test_shaft_warns_when_dynamic_efficiency_would_silently_cancel():
    # Shaft.efficiency scales the *net* shaft power, and a dynamic shaft's
    # steady equilibrium is exactly net == 0 -- so the factor cancels and
    # the mechanical loss disappears without any error. Belongs on the load
    # (ShaftLoad(efficiency=...)), where it scales a one-way power draw.
    comp, turb = _free_pair()
    with pytest.warns(UserWarning, match="efficiency scales the \\*net\\* shaft power"):
        Shaft(
            name="shaft", members=[comp, turb], signs=[-1.0, 1.0],
            efficiency=0.98, inertia=0.05, dynamic=True, N0=60000.0,
        )


def test_shaft_does_not_warn_for_dynamic_with_unity_efficiency():
    comp, turb = _free_pair()
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning becomes a failure
        Shaft(
            name="shaft", members=[comp, turb], signs=[-1.0, 1.0],
            efficiency=1.0, inertia=0.05, dynamic=True, N0=60000.0,
        )


def test_shaft_does_not_warn_for_static_efficiency():
    # In static mode net power is genuinely non-zero, so efficiency is real.
    comp, turb = _free_pair()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        Shaft(name="shaft", members=[comp, turb], signs=[-1.0, 1.0], efficiency=0.98)
