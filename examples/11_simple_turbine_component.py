"""SimpleTurbine: Source -> SimpleTurbine -> Sink, fixed expansion ratio +
isentropic efficiency (the turbine counterpart of SimpleCompressor).

Run: .venv/bin/python examples/11_simple_turbine_component.py
"""

from thermowave.components.simple_turbine import SimpleTurbine
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

src = Source(name="src", P=400000.0, T=1200.0, mdot=1.0)
turb = SimpleTurbine(name="t1", PR=2.5, eta_s=0.85, gamma=gamma)
snk = Sink(name="snk")

network = Network(fluid=air)
network.add_component(src)
network.add_component(turb)
network.add_component(snk)

network.connect(src, "out", turb, "in")
network.connect(turb, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=50, verbose=False)
result.print_report()
