"""Junction: Source -> Splitter -> {Valve A, Valve B} -> Merger -> Sink.

Demonstrates a genuinely branching network: the Network solver now carries a
per-node mass flow rate (not one global scalar), so a Splitter can send a
fixed fraction of the inlet flow down each branch and a Merger can sum two
independently-throttled branches back into one stream.

Run: .venv/bin/python examples/13_junction_splitter_merger.py
"""

from thermowave.components.junction import Junction
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.valve import Valve
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)

src = Source(name="src", P=200000.0, T=300.0, mdot=2.0)
splitter = Junction(name="split", n_inlets=1, n_outlets=2, split_fractions=[0.4, 0.6])
valve_a = Valve(name="va", D=0.1, K=3.0)
valve_b = Valve(name="vb", D=0.1, K=8.0)
merger = Junction(name="merge", n_inlets=2, n_outlets=1)
snk = Sink(name="snk")

network = Network(fluid=air)
for component in (src, splitter, valve_a, valve_b, merger, snk):
    network.add_component(component)

network.connect(src, "out", splitter, "in0")
network.connect(splitter, "out0", valve_a, "in")
network.connect(splitter, "out1", valve_b, "in")
network.connect(valve_a, "out", merger, "in0")
network.connect(valve_b, "out", merger, "in1")
network.connect(merger, "out0", snk, "in")

result = network.solve(tol=1e-8, max_iter=50, verbose=False)
result.print_report()
