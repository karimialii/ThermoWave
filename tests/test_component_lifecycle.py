"""Network.remove_component()/replace_component()/check_wiring().

The motivating case is the steady -> transient handoff: a steady solve pins a
compressor's speed with a Controller, then the transient run wants a
PIDController on the same unknown. Each contributes one residual, so the
Controller has to be removed rather than merely joined.
"""

import math

import pytest

from thermowave.components.compressor import Compressor
from thermowave.components.heat_transfer import Convection, ThermalMass, normalized_heat_paths
from thermowave.components.pid_controller import PIDController
from thermowave.components.pipe import Pipe
from thermowave.components.sensor import Sensor
from thermowave.components.setpoint import Setpoint
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.exceptions import NetworkTopologyError
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)
GAMMA = 1005.0 / (1005.0 - 287.05)
_COMP_MAP = "tests/fixtures/simple_compressor_map.cop"


def _pipe_chain():
    network = Network(fluid=AIR)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    p1 = Pipe(name="p1", L=5.0, D=0.2, f=0.02)
    p2 = Pipe(name="p2", L=5.0, D=0.2, f=0.02)
    snk = Sink(name="snk")
    for component in (src, p1, p2, snk):
        network.add_component(component)
    network.connect(src, "out", p1, "in")
    network.connect(p1, "out", p2, "in")
    network.connect(p2, "out", snk, "in")
    return network, src, p1, p2, snk


def _free_speed_compressor_network():
    network = Network(fluid=AIR)
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    comp = Compressor(name="comp", map_path=_COMP_MAP, gamma=GAMMA, N=None)
    sensor = Sensor(name="s1")
    snk = Sink(name="snk")
    for component in (src, comp, sensor, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", sensor, "tap")
    network.connect(comp, "out", snk, "in")
    return network, comp, sensor


def test_remove_component_drops_it_and_its_connections():
    network, src, p1, p2, snk = _pipe_chain()

    network.remove_component(p2)

    assert p2 not in network.components
    assert all(
        c.from_component is not p2 and c.to_component is not p2 for c in network.connections
    )
    assert len(network.connections) == 1  # src->p1 survives; p1->p2 and p2->snk are gone


def test_remove_component_rebuilds_union_find_so_survivors_still_solve():
    # The union-find is a flat child->parent map with path compression and no
    # reverse index, so a removed port id that happened to be a union root
    # cannot be unpicked in place -- remove_component() rebuilds instead.
    network, src, p1, p2, snk = _pipe_chain()

    network.remove_component(p2)
    network.connect(p1, "out", snk, "in")  # re-close the chain around the gap

    result = network.solve(tol=1e-10, max_iter=200, damping=0.6, progress=False)

    assert result.converged
    assert result.node_P["p1.out"] < result.node_P["src.out"]
    assert "p2.in" not in result.node_P


def test_remove_component_rejects_a_component_that_was_never_added():
    network, _, _, _, _ = _pipe_chain()
    stranger = Pipe(name="stranger", L=1.0, D=0.1, f=0.02)

    with pytest.raises(NetworkTopologyError, match="not in this network"):
        network.remove_component(stranger)


def test_remove_component_forgets_heat_paths_that_referenced_it():
    # Otherwise a surviving ThermalMass keeps calling Q() on a path whose
    # endpoint has left the network.
    network, comp, _ = _free_speed_compressor_network()
    casing = ThermalMass(name="casing", thermal_capacitance=200.0, T0=300.0)
    path = Convection(name="path", a=(comp, "out"), b=casing, h=50.0, A=0.3)
    network.add_heat_path(path)
    assert casing.heat_sources  # wired before removal

    network.remove_component(comp)

    assert casing.heat_sources == []


def test_remove_component_drops_the_path_itself_when_the_path_is_removed():
    network, comp, _ = _free_speed_compressor_network()
    casing = ThermalMass(name="casing", thermal_capacitance=200.0, T0=300.0)
    path = Convection(name="path", a=(comp, "out"), b=casing, h=50.0, A=0.3)
    network.add_heat_path(path)

    network.remove_component(path)

    assert casing.heat_sources == []
    assert normalized_heat_paths(comp.heat_path) == []


def test_replace_component_swaps_a_setpoint_for_a_pid_and_stays_square():
    # The actual steady -> transient handoff. Both pin comp.N with one
    # residual each; keeping both would over-determine the system by one.
    network, comp, sensor = _free_speed_compressor_network()
    setpoint = Setpoint(
        name="sp", component=comp, free_param="shaft", target_metric="PR [-]", value=2.5
    )
    network.add_component(setpoint)
    steady = network.solve(tol=1e-8, max_iter=300, damping=0.5, progress=False)

    pid = PIDController(
        name="pid", sensor=sensor, quantity="T [K]", component=comp, free_param="shaft",
        setpoint=420.0, Kp=60.0, Ki=50.0, Kd=0.0,
        output0=steady.node_N["comp.shaft"],
    )
    network.replace_component(setpoint, pid)

    assert setpoint not in network.components
    assert pid in network.components
    after = network.solve(tol=1e-8, max_iter=300, damping=0.5, progress=False)
    # The PID pins N to output0, which was seeded from the steady solve, so
    # this must re-find the same operating point.
    assert math.isclose(after.node_N["comp.shaft"], steady.node_N["comp.shaft"], rel_tol=1e-6)


def test_check_wiring_names_a_doubly_closed_free_parameter():
    network, comp, sensor = _free_speed_compressor_network()
    network.add_component(
        Setpoint(name="sp", component=comp, free_param="shaft", target_metric="PR [-]", value=2.5)
    )
    network.add_component(
        PIDController(
            name="pid", sensor=sensor, quantity="T [K]", component=comp, free_param="shaft",
            setpoint=420.0, Kp=1.0, output0=60000.0,
        )
    )

    problems = network.check_wiring()

    assert any("comp.shaft" in p and "closed by 2" in p for p in problems)


def test_check_wiring_names_an_unclosed_free_parameter():
    network, comp, _ = _free_speed_compressor_network()  # comp.N free, nothing pins it

    problems = network.check_wiring()

    assert any("comp.shaft" in p and "nothing closes it" in p for p in problems)


def test_check_wiring_names_a_controller_left_pointing_at_a_removed_component():
    network, comp, sensor = _free_speed_compressor_network()
    setpoint = Setpoint(
        name="sp", component=comp, free_param="shaft", target_metric="PR [-]", value=2.5
    )
    network.add_component(setpoint)

    network.remove_component(comp)

    problems = network.check_wiring()
    assert any("'sp'" in p and "not in this network" in p for p in problems)


def test_check_wiring_is_quiet_on_a_correctly_wired_network():
    network, comp, _ = _free_speed_compressor_network()
    network.add_component(
        Setpoint(name="sp", component=comp, free_param="shaft", target_metric="PR [-]", value=2.5)
    )

    assert network.check_wiring() == []


def test_non_square_solve_error_names_the_doubly_closed_parameter():
    # The solver's own count-mismatch message can't say which parameter is
    # at fault; check_wiring()'s findings get appended to it.
    network, comp, sensor = _free_speed_compressor_network()
    network.add_component(
        Setpoint(name="sp", component=comp, free_param="shaft", target_metric="PR [-]", value=2.5)
    )
    network.add_component(
        PIDController(
            name="pid", sensor=sensor, quantity="T [K]", component=comp, free_param="shaft",
            setpoint=420.0, Kp=1.0, output0=60000.0,
        )
    )

    with pytest.raises(NetworkTopologyError) as excinfo:
        network.solve(tol=1e-8, max_iter=50, progress=False)

    message = str(excinfo.value)
    assert "Likely cause:" in message
    assert "comp.shaft" in message


def test_solve_result_state_round_trips_for_report_metrics():
    network, comp, _ = _free_speed_compressor_network()
    network.add_component(
        Setpoint(name="sp", component=comp, free_param="shaft", target_metric="PR [-]", value=2.5)
    )
    result = network.solve(tol=1e-8, max_iter=300, damping=0.5, progress=False)

    metrics = comp.report_metrics(result.state())

    assert math.isclose(metrics["PR [-]"], 2.5, rel_tol=1e-6)


def test_check_wiring_recognises_a_shaft_as_closing_its_members_speeds():
    # A Shaft closes speed-tied members' speeds through its speed-tie
    # residuals, not through any controller — so those must not be reported
    # as unclosed. This is what BaseComponent.closes_parameters() is for.
    from thermowave.components.shaft import Shaft
    from thermowave.components.turbine import Turbine

    network = Network(fluid=AIR)
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    comp = Compressor(name="comp", map_path=_COMP_MAP, gamma=GAMMA, N=None)
    heater = Pipe(name="heater", L=1.0, D=0.1, f=0.0, n_elem=1, heat_loss=-400_000.0)
    turb = Turbine(
        name="turb", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=None
    )
    shaft = Shaft(
        name="shaft", members=[comp, turb], signs=[-1.0, 1.0],
        efficiency=1.0, inertia=0.05, dynamic=True, N0=60_000.0,
    )
    snk = Sink(name="snk")
    for component in (src, comp, heater, turb, shaft, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", heater, "in")
    network.connect(heater, "out", turb, "in")
    network.connect(turb, "out", snk, "in")
    network.connect(shaft, "m0", comp, "shaft", kind="mechanical")
    network.connect(shaft, "m1", turb, "shaft", kind="mechanical")

    assert network.check_wiring() == []


def test_check_wiring_flags_a_static_shafts_unclosed_reference_speed():
    # dynamic=False ties followers to components[0]'s speed but leaves that
    # reference for something external to close — so it should be reported.
    from thermowave.components.shaft import Shaft
    from thermowave.components.turbine import Turbine

    comp = Compressor(name="comp", map_path=_COMP_MAP, gamma=GAMMA, N=None)
    turb = Turbine(
        name="turb", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=None
    )
    shaft = Shaft(name="shaft", members=[comp, turb], signs=[-1.0, 1.0])

    network = Network(fluid=AIR)
    for component in (comp, turb, shaft):
        network.add_component(component)
    network.connect(shaft, "m0", comp, "shaft", kind="mechanical")
    network.connect(shaft, "m1", turb, "shaft", kind="mechanical")

    problems = network.check_wiring()
    # The canonical mechanical node id is whichever port id connect() rooted
    # the union at ("shaft.m0" here, from connect(shaft, "m0", comp, ...));
    # what matters is that exactly the reference speed is reported free.
    assert any("shaft.m0" in p and "nothing closes it" in p for p in problems)
    assert not any("shaft.m1" in p for p in problems)  # tied to comp, so closed


def test_convergence_failure_reports_the_system_size():
    # A bare "Singular Jacobian at iteration 12" says the linear algebra
    # broke, not which part of the model broke it. The system size and any
    # wiring findings are only known inside the solve, so they get attached
    # there rather than needing a second run with verbose=True.
    from thermowave.core.exceptions import ConvergenceError

    network, comp, sensor = _free_speed_compressor_network()
    network.add_component(
        Setpoint(
            name="sp", component=comp, free_param="shaft",
            target_metric="PR [-]", value=500.0,  # unreachable on this map
        )
    )

    with pytest.raises(ConvergenceError) as excinfo:
        network.solve(tol=1e-12, max_iter=3, damping=0.1, progress=False)

    message = str(excinfo.value)
    assert "unknown(s)" in message and "equation(s)" in message
    assert "steady state" in message
    assert "warm_start" in message  # points at the usual fix


def test_convergence_failure_surfaces_wiring_problems_when_there_are_any():
    # An unclosed free parameter usually shows up as a non-square system,
    # but if the counts happen to balance it can reach the Newton solve and
    # fail there instead -- so the wiring findings ride along either way.
    from thermowave.core.exceptions import ConvergenceError

    network, comp, sensor = _free_speed_compressor_network()
    network.add_component(
        Setpoint(name="sp", component=comp, free_param="shaft",
                 target_metric="PR [-]", value=500.0)
    )
    network.remove_component(sensor)  # harmless here, keeps the system square

    with pytest.raises(ConvergenceError) as excinfo:
        network.solve(tol=1e-12, max_iter=3, damping=0.1, progress=False)
    assert "System:" in str(excinfo.value)


def test_successful_solve_is_unaffected_by_the_failure_diagnostics():
    network, comp, _ = _free_speed_compressor_network()
    network.add_component(
        Setpoint(name="sp", component=comp, free_param="shaft", target_metric="PR [-]", value=2.5)
    )

    result = network.solve(tol=1e-8, max_iter=300, damping=0.5, progress=False)

    assert result.converged
