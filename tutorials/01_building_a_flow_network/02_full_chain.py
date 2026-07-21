"""Tutorial 1, step 2: a fuller chain — Source -> Compressor -> Pipe -> Valve -> Sink.

Step 1 built the two-node skeleton (Source -> Sink). Here every component in
between contributes its own physics to the same solve: SimpleCompressor
raises pressure via a fixed ratio and isentropic efficiency, Pipe adds
friction pressure drop and heat loss over its length, and Valve adds a
K-factor pressure drop scaled by how open it is. None of them needed any new
network-building mechanics beyond what step 1 already showed — add_component()
+ connect() is the whole pattern, repeated once per component.

What you'll learn:
  - Chaining four components is just four connect() calls; Network.solve()
    doesn't care how long the chain is.
  - Every component reads the same shared fluid (air) via the network's
    (P, h) state — SimpleCompressor and Pipe and Valve never see each other
    directly, only the node state between them.
  - result.print_report() groups components by category (TurboMachinery for
    the compressor) and always prints every node's state at the end.

Run: .venv/bin/python tutorials/01_building_a_flow_network/02_full_chain.py
"""

from thermowave.components.pipe import Pipe
from thermowave.components.simple_compressor import SimpleCompressor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.valve import Valve
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
comp = SimpleCompressor(name="c1", PR=2.5, eta_s=0.8, gamma=gamma)
pipe = Pipe(name="p1", L=6.0, D=0.2, f=0.02, n_elem=5, heat_loss=3000.0)
valve = Valve(name="v1", D=0.2, K=2.0, opening=0.7)
snk = Sink(name="snk")

network = Network(fluid=air)
for component in (src, comp, pipe, valve, snk):
    network.add_component(component)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", pipe, "in")
network.connect(pipe, "out", valve, "in")
network.connect(valve, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=50, verbose=True)
result.print_report()

print(
    "\nThis is as far as a pure flow chain goes. Building an actual gas "
    "turbine cycle — with a compressor and turbine sharing one shaft, a "
    "combustor, and eventually a recuperator — is tutorial 2: "
    "tutorials/02_gas_turbine_cycle/."
)
