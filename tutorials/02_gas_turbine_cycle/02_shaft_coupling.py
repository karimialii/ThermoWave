"""Tutorial 2, step 2: putting the compressor and turbine on one shaft.

Step 1 gave the compressor and turbine independent, fixed speeds — a
convenient fiction, but not how a real single-shaft engine works: there's
one physical rotor, so the compressor and turbine turn at exactly the same
speed, and whatever speed that settles at is determined by the whole
cycle's torque balance, not chosen in advance for each machine separately.

Shaft is how ThermoWave expresses that constraint. Both machines are left
with N=None (free, solved-for unknowns) instead of a fixed number, and a
Shaft ties them to the same speed (signs=[-1.0, 1.0] tells it the compressor
draws power from the shaft while the turbine delivers power to it). That
introduces two free unknowns (comp.N, turb.N) but only one physical degree
of freedom (they're tied together) — a Setpoint on the compressor's own N,
still pinning it to the same 70000 rev/min step 1 used directly, supplies
the one closing equation.

What you'll learn:
  - Shaft(components=[...], signs=[...]) mechanically couples any number of
    rotating components' speeds; efficiency scales for bearing/windage
    losses on top of that coupling.
  - Setpoint(component, free_param, target_metric, value) is the general
    mechanism for pinning any component's free parameter to a target metric
    — not specific to shafts. It's what closes the one remaining degree of
    freedom Shaft's coupling leaves open.
  - Compare this run's TurboMachinery table to step 1's: same physics, same
    N, but now shaft.power [W] reports the actual net mechanical
    surplus (turbine output minus compressor draw, minus bearing losses)
    instead of you having to subtract the two numbers yourself.

Run: .venv/bin/python tutorials/02_gas_turbine_cycle/02_shaft_coupling.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.setpoint import Setpoint
from thermowave.components.shaft import Shaft
from thermowave.components.simple_combustor import SimpleCombustor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

target_N = 70000.0  # rev/min — same value step 1 fixed directly

src = Source(name="src", P=101325.0, T=288.15, mdot=0.8)
comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=gamma, N=None)
combustor = SimpleCombustor(name="cc", LHV=50.0e6, PR=0.96, efficiency=0.99, mdot_fuel=0.012)
turb = Turbine(name="turb", map_path="T100 Turb.tur", gamma=gamma, N=None)
shaft = Shaft(name="shaft", components=[comp, turb], signs=[-1.0, 1.0], efficiency=0.98)
setpoint = Setpoint(
    name="sp_N", component=comp, free_param="N", target_metric="N [rev/min]", value=target_N
)
snk = Sink(name="snk")

network = Network(fluid=air)
for component in (src, comp, combustor, turb, shaft, setpoint, snk):
    network.add_component(component)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", combustor, "in")
network.connect(combustor, "out", turb, "in")
network.connect(turb, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=100, damping=0.5, verbose=True)
result.print_report()

print(
    "\nBoth machines converged to the same N the shaft coupling now enforces "
    "— compare shaft.power [W] to step 1's turb-minus-comp power by hand.\n"
    "Step 3 adds a recuperator (preheating compressor discharge with turbine "
    "exhaust heat) and a physical fuel supply line, completing the cycle."
)
