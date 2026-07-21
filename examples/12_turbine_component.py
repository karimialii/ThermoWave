"""Turbine: Source -> Turbine -> Sink, driven by a Flownex-style
characteristic map (T100 Turb.tur — a real Turbec T100 microturbine turbine
map, iso-speed curves of expansion ratio and efficiency vs non-dimensional
mass flow).

Unlike SimpleTurbine (fixed PR + eta_s), Turbine looks PR and eta_s up
from the map at the current corrected speed (derived from N and the current
inlet temperature) and corrected mass flow each solver iteration.

Run: .venv/bin/python examples/12_turbine_component.py
"""

from thermowave.components.turbine import Turbine
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

N = 70000.0  # real shaft speed [rev/min]; Turbine derives corrected speed itself

src = Source(name="src", P=400000.0, T=1150.0, mdot=0.8)
turb = Turbine(name="t1", map_path="T100 Turb.tur", N=N, gamma=gamma)
snk = Sink(name="snk")

network = Network(fluid=air)
network.add_component(src)
network.add_component(turb)
network.add_component(snk)

network.connect(src, "out", turb, "in")
network.connect(turb, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=50, verbose=True)
result.print_report()
