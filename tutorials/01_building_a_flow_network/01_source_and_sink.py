"""Tutorial 1, step 1: the smallest possible network.

Every ThermoWave network needs three things: a fluid model, at least one
component that fixes a boundary condition (a Source pins pressure,
temperature, and mass flow), and a component that terminates the flow (a
Sink). This is the minimum viable network — two components, one connection.

What you'll learn:
  - Network(fluid=...) owns the fluid model every component in it shares.
  - add_component() registers a component; connect() wires one component's
    named port to another's.
  - network.solve() runs the Newton solve; result.print_report() shows the
    solved (P, T, h, mdot) at every node.

Since the Sink has no outlet, there's nothing to compute here — the network
just reports back exactly what the Source fixed. The point is the mechanics
of assembling a network, not the physics yet; that starts in step 2.

Run: .venv/bin/python tutorials/01_building_a_flow_network/01_source_and_sink.py
"""

from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)

src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
snk = Sink(name="snk")

network = Network(fluid=air)
network.add_component(src)
network.add_component(snk)

network.connect(src, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=50, verbose=True)
result.print_report()

print(
    "\nNext: 02_full_chain.py adds a compressor, a pipe, and a valve between "
    "the Source and Sink, so the solve actually has physics to do."
)
