"""Tutorial 2, step 1: a basic Brayton cycle — Source -> Compressor -> combustor -> Turbine -> Sink.

The core of every gas turbine: air is compressed, heated by burning fuel,
and expanded through a turbine that extracts more work than the compressor
consumed (that surplus is the machine's useful output). Here the compressor
and turbine speeds are both given directly — as if an external generator
holds shaft speed fixed — so there's no coupling between them yet: each one
independently follows its own performance map at that speed. Step 2 removes
that assumption and actually ties them together on one physical shaft.

Compressor and Turbine (not SimpleCompressor/SimpleTurbine) are used here:
instead of a fixed pressure ratio and efficiency, they look up pressure
ratio and efficiency from a real compressor/turbine performance map
(T100 Comp.cop / T100 Turb.tur, corrected-mass-flow vs. corrected-speed
maps for a small commercial microturbine) given the current corrected speed
and mass flow. SimpleCombustor adds fuel mass and heat via a fixed lower
heating value (LHV), rather than genuine combustion chemistry.

What you'll learn:
  - Map-based Compressor/Turbine take N [rev/min] directly, the same way
    SimpleCompressor takes a fixed PR.
  - SimpleCombustor's mdot_fuel can be a plain fixed input (as here) or left
    free — that's step 3's territory.
  - gamma (needed by the maps' isentropic relations) can be derived from any
    BaseFluid via BaseFluid.gamma(), so it doesn't have to be hardcoded per
    fluid.

Run: .venv/bin/python tutorials/02_gas_turbine_cycle/01_basic_brayton_cycle.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.simple_combustor import SimpleCombustor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

shaft_N = 70000.0  # rev/min, held fixed for both machines (no shaft coupling yet)

src = Source(name="src", P=101325.0, T=288.15, mdot=0.8)
comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=gamma, N=shaft_N)
combustor = SimpleCombustor(name="cc", LHV=50.0e6, PR=0.96, efficiency=0.99, mdot_fuel=0.012)
turb = Turbine(name="turb", map_path="T100 Turb.tur", gamma=gamma, N=shaft_N)
snk = Sink(name="snk")

network = Network(fluid=air)
for component in (src, comp, combustor, turb, snk):
    network.add_component(component)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", combustor, "in")
network.connect(combustor, "out", turb, "in")
network.connect(turb, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=100, damping=0.5, verbose=True)
result.print_report()

print(
    "\nNotice compressor and turbine power in the report aren't equal —\n"
    "nothing here ties them together yet. Step 2 puts them on one shaft."
)
