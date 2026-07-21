"""Tutorial 3, step 2: watching the shaft actually spool up over time.

Step 1's steady solve found the rotor's equilibrium speed but can't show
*how it gets there* — a steady solve of a differential-state system only
ever sees the state where nothing is changing. Network.solve_transient()
is the tool for genuine time-domain behavior: it integrates every
differential state (here, just shaft.N) forward via backward-Euler,
re-solving the rest of the network's algebraic state at each timestep.

This starts the shaft 10000 rev/min below step 1's equilibrium (net
positive torque — the turbine is out-producing what the compressor draws at
that speed) and watches it spool up. The `initial=` SolveResult handed to
solve_transient() is a deliberately off-equilibrium t=0 condition — pass
nothing and solve_transient() runs an ordinary Network.solve() first, which
(since N closes via "derivative == 0") is already sitting exactly at
equilibrium, so nothing would move.

What you'll learn:
  - solve_transient(duration, dt, initial=...) returns a TransientResult:
    one SolveResult per timestep (history.steps), plus a
    name -> value-list history for every differential state
    (history.diff_history["shaft.N"]).
  - Building an off-equilibrium initial condition: copy a converged
    SolveResult and edit its .params dict directly.
  - history.plot(...) renders any component's report_metrics() over time —
    the same call signature works for any component, not just Shaft.

Run: .venv/bin/python tutorials/03_transient_rotor_dynamics/02_spool_up_transient.py
"""

import copy

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

# Steady-state equilibrium speed (same as step 1), then a deliberately
# off-equilibrium t=0 condition 10000 rev/min below it.
equilibrium = network.solve(tol=1e-8, max_iter=400, damping=0.3)
print(f"Equilibrium speed: {equilibrium.params['shaft.N']:.1f} rev/min")

initial = copy.copy(equilibrium)
initial.params = dict(equilibrium.params)
initial.params["shaft.N"] = equilibrium.params["shaft.N"] - 10000.0

history = network.solve_transient(
    duration=1.0, dt=0.05, initial=initial, tol=1e-8, max_iter=400, damping=0.3,
)

print(f"\n{'t [s]':>8}  {'N [rev/min]':>12}")
for t, N in zip(history.times, history.diff_history["shaft.N"]):
    print(f"{t:8.2f}  {N:12.1f}")

print(f"\nFinal state at t = {history.times[-1]:.2f} s:")
history.steps[-1].print_report()

# Plotting requires the 'plot' extra: pip install thermowave[plot]. show=False
# + save_path here so this example stays non-interactive when run headless/in
# CI — drop both (or just show=True) to pop up a window interactively instead.
history.plot(
    (shaft, "N [rev/min]"), ylabel="N [rev/min]", title="Shaft spool-up",
    show=False, save_path="spool_up_speed.png",
)
history.plot(
    (comp, "power [W]"), (turb, "power [W]"),
    ylabel="power [W]", title="Compressor vs. turbine power",
    show=False, save_path="spool_up_power.png",
)
print("\nSaved spool_up_speed.png and spool_up_power.png")

print(
    "\nNext: tutorials/04_closed_loop_control/ builds real PID control loops "
    "on top of a steady-state cycle — a different kind of time-varying "
    "behavior than rotor dynamics: the setpoint changes, not the physics."
)
