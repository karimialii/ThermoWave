"""Compressor: Source -> Compressor -> Sink, driven by a Flownex-style
characteristic map (T100 Comp.cop — a real Turbec T100 microturbine
compressor map, iso-speed curves of pressure ratio and efficiency vs
non-dimensional mass flow).

Unlike SimpleCompressor (fixed PR + eta_s), Compressor looks PR and eta_s
up from the map at the current corrected speed (from N) and corrected mass
flow (from the current inlet state) each solver iteration.

Run: .venv/bin/python examples/10_compressor_component.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

T_in = 288.15
N = 70000.0  # real shaft speed [rev/min]; Compressor derives corrected speed itself

src = Source(name="src", P=101325.0, T=T_in, mdot=0.8)
comp = Compressor(name="c1", map_path="T100 Comp.cop", N=N, gamma=gamma)
snk = Sink(name="snk")

network = Network(fluid=air)
network.add_component(src)
network.add_component(comp)
network.add_component(snk)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=50, verbose=False)
result.print_report()
