"""Tutorial 3, step 1: a shaft that finds its own equilibrium speed.

Tutorial 2's shaft coupling used Shaft's default steady-state mode: N is a
free algebraic unknown, closed by an external Setpoint pinning it to a
chosen value. That models a control loop holding speed, but says nothing
about the rotor's own dynamics — the shaft has no mass, no inertia, no
sense of "spinning up" or "slowing down" over time.

Shaft(dynamic=True, inertia=...) changes what closes N: instead of a
Setpoint, N becomes a genuine differential state (see
BaseComponent.differential_parameters()), driven by
    d(N)/dt = (net_shaft_power / omega) / inertia
A plain Network.solve() (steady state) then finds N by requiring this
derivative to be exactly zero — i.e. the speed at which the compressor's
power draw and the turbine's power output are in torque balance, the
rotor's genuine equilibrium, discovered by the physics rather than chosen
by a Setpoint. (The heater here is a Pipe with a fixed negative heat_loss —
a massless heat-addition placeholder, not a real combustor — so this step
stays focused on the shaft; tutorial 2's SimpleCombustor cycle returns in
step 2's own combustor-based variant if you want the two combined.)

What you'll learn:
  - dynamic=True + inertia turns Shaft's speed into differential state,
    removing the need for any Setpoint/Controller to close it.
  - A steady Network.solve() of a differential-state system finds the
    equilibrium (derivative == 0), not a snapshot in time — genuine
    time-domain behavior needs solve_transient(), which step 2 covers.

Run: .venv/bin/python tutorials/03_transient_rotor_dynamics/01_dynamic_shaft_equilibrium.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.pipe import Pipe
from thermowave.components.shaft import Shaft
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=gamma, N=None)
heater = Pipe(name="heater", L=1.0, D=0.1, f=0.0, n_elem=1, heat_loss=-431000.0)
turb = Turbine(name="turb", map_path="T100 Turb.tur", gamma=gamma, N=None)
shaft = Shaft(
    name="shaft", components=[comp, turb], signs=[-1.0, 1.0], efficiency=0.98,
    inertia=0.05, dynamic=True, N0=55000.0,
)
snk = Sink(name="snk")

network = Network(fluid=air)
for component in (src, comp, heater, turb, shaft, snk):
    network.add_component(component)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", heater, "in")
network.connect(heater, "out", turb, "in")
network.connect(turb, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=400, damping=0.3, verbose=True)
result.print_report()

print(f"\nEquilibrium speed found by the solver: {result.params['shaft.N']:.1f} rev/min")
print(
    "\nNo Setpoint anywhere in this network, yet the shaft still landed on a "
    "specific speed — that's the torque-balance equilibrium, not a chosen "
    "value. Step 2 starts below this equilibrium and watches the shaft spool "
    "up to it over time."
)
