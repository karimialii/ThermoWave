import math
from pathlib import Path

import pytest

from thermowave.components.compressor import Compressor
from thermowave.components.schedule import Schedule
from thermowave.components.setpoint import Setpoint
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

_MAP_PATH = str(Path(__file__).parent / "fixtures" / "simple_compressor_map.cop")

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _Target:
    def __init__(self, value):
        self.value = value


def test_schedule_rejects_fewer_than_two_breakpoints():
    target = _Target(0.0)
    with pytest.raises(ValueError, match="at least two"):
        Schedule(name="sch1", target=target, attr="value", breakpoints=[(0.0, 1.0)])


def test_schedule_rejects_non_increasing_times():
    target = _Target(0.0)
    with pytest.raises(ValueError, match="strictly increasing"):
        Schedule(
            name="sch1", target=target, attr="value",
            breakpoints=[(0.0, 1.0), (5.0, 2.0), (3.0, 3.0)],
        )


def test_schedule_rejects_duplicate_times():
    target = _Target(0.0)
    with pytest.raises(ValueError, match="strictly increasing"):
        Schedule(
            name="sch1", target=target, attr="value",
            breakpoints=[(0.0, 1.0), (0.0, 2.0)],
        )


def test_schedule_rejects_unknown_interpolation():
    target = _Target(0.0)
    with pytest.raises(ValueError, match="interpolation"):
        Schedule(
            name="sch1", target=target, attr="value",
            breakpoints=[(0.0, 1.0), (1.0, 2.0)], interpolation="cubic",
        )


def test_schedule_rejects_target_missing_attr():
    target = object()
    with pytest.raises(ValueError, match="doesn't exist"):
        Schedule(name="sch1", target=target, attr="value", breakpoints=[(0.0, 1.0), (1.0, 2.0)])


def test_schedule_seeds_target_attr_at_construction_to_t0_value():
    target = _Target(0.0)
    Schedule(name="sch1", target=target, attr="value", breakpoints=[(0.0, 5.0), (10.0, 15.0)])
    assert target.value == 5.0


def test_schedule_ports_and_residuals_are_empty():
    target = _Target(0.0)
    sch = Schedule(name="sch1", target=target, attr="value", breakpoints=[(0.0, 1.0), (1.0, 2.0)])
    assert sch.ports() == {}
    assert sch.residuals(state=None) == []


def test_schedule_linear_interpolation_between_breakpoints():
    target = _Target(0.0)
    sch = Schedule(
        name="sch1", target=target, attr="value",
        breakpoints=[(0.0, 0.0), (10.0, 100.0)],
    )
    sch.step(state=None, dt=2.5)
    assert math.isclose(target.value, 25.0)
    sch.step(state=None, dt=2.5)
    assert math.isclose(target.value, 50.0)


def test_schedule_step_interpolation_holds_earlier_value_until_next_breakpoint():
    target = _Target(0.0)
    sch = Schedule(
        name="sch1", target=target, attr="value",
        breakpoints=[(0.0, 1.0), (10.0, 2.0), (20.0, 3.0)],
        interpolation="step",
    )
    sch.step(state=None, dt=5.0)
    assert target.value == 1.0
    sch.step(state=None, dt=5.0)  # now at t=10.0, exactly on the second breakpoint
    assert target.value == 2.0
    sch.step(state=None, dt=5.0)
    assert target.value == 2.0


def test_schedule_holds_flat_before_first_and_after_last_breakpoint():
    target = _Target(0.0)
    sch = Schedule(
        name="sch1", target=target, attr="value",
        breakpoints=[(5.0, 1.0), (10.0, 2.0)],
    )
    assert target.value == 1.0  # t=0 is before the first breakpoint
    sch.step(state=None, dt=100.0)  # far past the last breakpoint
    assert target.value == 2.0


def test_schedule_report_metrics_reflects_current_target_value():
    target = _Target(0.0)
    sch = Schedule(name="sch1", target=target, attr="value", breakpoints=[(0.0, 3.0), (1.0, 3.0)])
    assert sch.report_metrics(state=None) == {"target [-]": 3.0}


def test_schedule_drives_compressor_pr_setpoint_over_transient():
    gamma = 1005.0 / (1005.0 - 287.05)
    src = Source(name="src", P=101325.0, T=300.0, mdot=2.0 * 1.01325 / math.sqrt(300.0))
    comp = Compressor(name="c1", map_path=_MAP_PATH, gamma=gamma)  # N left free
    sp = Setpoint(name="sp1", component=comp, free_param="N", target_metric="PR [-]", value=3.0)
    sch = Schedule(
        name="sch1", target=sp, attr="value",
        breakpoints=[(0.0, 3.0), (5.0, 3.5)],
    )
    snk = Sink(name="snk")

    network = Network(fluid=AIR)
    for component in (src, comp, sp, sch, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", snk, "in")

    # step() only updates the schedule's target after each solve (same lag as
    # PIDController.step() — see Schedule.step()'s docstring), so run one
    # extra step past the last breakpoint to let the final value land.
    history = network.solve_transient(duration=6.0, dt=1.0, tol=1e-8, max_iter=100)

    def _pr(step):
        state = NetworkState(
            fluid=step.fluid, node_P=step.node_P, node_h=step.node_h,
            node_mdot=step.node_mdot, params=step.params,
        )
        return comp.report_metrics(state)["PR [-]"]

    assert math.isclose(_pr(history.steps[0]), 3.0, rel_tol=1e-4)
    assert math.isclose(_pr(history.steps[-1]), 3.5, rel_tol=1e-4)
