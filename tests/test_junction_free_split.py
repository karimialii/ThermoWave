"""End-to-end: Junction's free split fractions + Controller's live (callable)
target, together closing a genuinely implicit flow split -- the pattern a
looped/recompression network needs (see Junction's own docstring) that a
fixed split_fractions input can't represent.

The network: one inlet splits into two pipes of different resistance, then
re-merges. Real flow physics don't pin a unique split on their own here (the
merge just mixes whatever arrives, no equal-pressure constraint of its own
-- see Junction's docstring), so this drives the split with an explicit,
physically meaningful design rule instead: the two branches must reach the
SAME pressure before merging (the real constraint an actual looped pipe
network enforces, and exactly the "unknown, pressure-balance-determined
split" case Junction's fixed split_fractions can't model on its own). A
Controller closes the free fraction by reading one branch's pressure via a
Sensor and matching it live against the other branch's Sensor reading --
the callable `value` this test exercises.
"""

import math

from thermowave.components.controller import Controller
from thermowave.components.junction import Junction
from thermowave.components.pipe import Pipe
from thermowave.components.sensor import Sensor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def _build_network(split_frac0):
    """split_frac0: fixed float (baseline network) or None (free, driven by
    a Controller matching both branches' pressure)."""
    net = Network(fluid=AIR)

    src = Source(name="src", P=300000.0, T=300.0, mdot=0.4)
    splitter = Junction(
        name="splitter", n_inlets=1, n_outlets=2, split_fractions=[split_frac0, None]
    )
    # Branch A: narrow pipe, more resistance per unit mass flow.
    pipe_a = Pipe(name="pipe_a", L=50.0, D=0.08, roughness=4.5e-5, mu=1.8e-5)
    sensor_a = Sensor(name="sensor_a")
    # Branch B: wide pipe, less resistance.
    pipe_b = Pipe(name="pipe_b", L=50.0, D=0.12, roughness=4.5e-5, mu=1.8e-5)
    sensor_b = Sensor(name="sensor_b")
    merger = Junction(name="merger", n_inlets=2, n_outlets=1)
    snk = Sink(name="snk")

    for c in (src, splitter, pipe_a, sensor_a, pipe_b, sensor_b, merger, snk):
        net.add_component(c)
    net.connect(src, "out", splitter, "in0")
    net.connect(splitter, "out0", pipe_a, "in")
    net.connect(pipe_a, "out", sensor_a, "tap")
    net.connect(sensor_a, "tap", merger, "in0")
    net.connect(splitter, "out1", pipe_b, "in")
    net.connect(pipe_b, "out", sensor_b, "tap")
    net.connect(sensor_b, "tap", merger, "in1")
    net.connect(merger, "out0", snk, "in")

    return net, splitter, sensor_a, sensor_b


def test_free_split_fraction_declared_and_solved_for():
    net, splitter, sensor_a, sensor_b = _build_network(split_frac0=None)
    ctrl = Controller(
        name="match_p", sensor=sensor_a, quantity="P [Pa]",
        component=splitter, free_param="frac0",
        value=lambda s: sensor_b.report_metrics(s)["P [Pa]"],
    )
    net.add_component(ctrl)

    result = net.solve(tol=1e-8, max_iter=200, progress=False)
    assert result.converged

    state = result.state()
    P_a = sensor_a.report_metrics(state)["P [Pa]"]
    P_b = sensor_b.report_metrics(state)["P [Pa]"]
    assert math.isclose(P_a, P_b, rel_tol=1e-6)

    # The narrower, higher-resistance branch (A) should carry LESS mass
    # flow than the wider one (B) at the pressure-matched solution.
    frac0 = state.param("splitter.frac0")
    assert 0.0 < frac0 < 0.5


def test_fixed_symmetric_split_does_not_match_branch_pressures():
    # Sanity check on the test setup itself: an even (fixed) 50/50 split
    # between two DIFFERENT-resistance pipes should NOT land on equal
    # pressures -- confirming the free-split solve above is doing real work,
    # not just reproducing what a fixed split already gives for free.
    net, splitter, sensor_a, sensor_b = _build_network(split_frac0=0.5)
    result = net.solve(tol=1e-8, max_iter=200, progress=False)
    assert result.converged

    state = result.state()
    P_a = sensor_a.report_metrics(state)["P [Pa]"]
    P_b = sensor_b.report_metrics(state)["P [Pa]"]
    assert not math.isclose(P_a, P_b, rel_tol=1e-3)
