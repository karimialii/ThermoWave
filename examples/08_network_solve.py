"""Network + Solver: assemble Source -> Pipe -> Valve -> Sink and solve end-to-end.

Components are added with add_component(), then wired together with
connect(from_component, from_port, to_component, to_port) — a typed
connection (kind="flow" for now; mechanical/signal/heat kinds are planned).

Demonstrates the verbose Newton-Raphson iteration reporter and the
SolveResult.print_report() table.

Run: .venv/bin/python examples/08_network_solve.py
"""

from thermowave.components.pipe import Pipe
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.valve import Valve
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)

src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
pipe = Pipe(name="p1", L=6.0, D=0.2, f=0.02, n_elem=10, heat_loss=0)
valve = Valve(name="v1", D=0.2, K=2.0, opening=.1)

snk = Sink(name="snk")

network = Network(fluid=air)
network.add_component(src)
network.add_component(pipe)
network.add_component(valve)
network.add_component(snk)

network.connect(src, "out", pipe, "in")
network.connect(pipe, "out", valve, "in")
network.connect(valve, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=50, verbose=False)
print()
result.print_report()
