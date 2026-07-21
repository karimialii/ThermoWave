"""Tutorial 5, step 3: closing the loop with a condenser and computing cycle efficiency.

A SimpleCondenser now takes the turbine's wet exhaust from step 2 and
rejects heat to bring it back to saturated liquid — the same outlet-state
targeting SimpleEvaporator uses (0.0 quality here instead of superheat),
just running in reverse. Its outlet state matches step 1's feedwater Source
almost exactly (both are liquid water at condenser pressure), which is what
makes this network the *unrolled* equivalent of a real closed Rankine loop:
a Source pins the feedwater state and a Sink terminates the exhaust (the
solver needs a fixed boundary node, and a genuinely closed recycle loop is
out of scope — see the main README's fluid-propagation section), but
reading the condenser outlet against the feed inlet shows they're
essentially the same state, so the four corner enthalpies below give the
real cycle numbers.

What you'll learn:
  - SimpleCondenser(outlet_quality=0.0) mirrors SimpleEvaporator, just
    rejecting heat instead of adding it.
  - Cycle efficiency from four corner enthalpies: net work (turbine output
    minus pump input) divided by heat input — the standard Rankine-cycle
    calculation, done here directly from the solved node enthalpies rather
    than a dedicated "cycle efficiency" method, since it's just arithmetic
    once every state point is known.

Requires the 'coolprop' extra: pip install thermowave[coolprop]

Run: .venv/bin/python tutorials/05_rankine_steam_cycle/03_full_cycle_with_condenser.py
"""

from thermowave.components.pump import Pump
from thermowave.components.simple_condenser import SimpleCondenser
from thermowave.components.simple_evaporator import SimpleEvaporator
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.steam_turbine import SteamTurbine
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.real_fluid import CoolPropFluid

water = CoolPropFluid(name="Water")

P_cond = 1.0e4    # 0.1 bar condenser pressure
P_boiler = 4.0e6  # 40 bar boiler pressure
mdot = 10.0       # kg/s

feed = Source(name="feed", P=P_cond, T=315.0, mdot=mdot)
pump = Pump(name="pump", P_out=P_boiler, eta=0.75)
boiler = SimpleEvaporator(name="boiler", superheat=150.0)
turbine = SteamTurbine(name="turbine", P_out=P_cond, eta_s=0.85)
condenser = SimpleCondenser(name="condenser", outlet_quality=0.0)
exhaust = Sink(name="exhaust")

network = Network(fluid=water)
for component in (feed, pump, boiler, turbine, condenser, exhaust):
    network.add_component(component)

network.connect(feed, "out", pump, "in")
network.connect(pump, "out", boiler, "in")
network.connect(boiler, "out", turbine, "in")
network.connect(turbine, "out", condenser, "in")
network.connect(condenser, "out", exhaust, "in")

result = network.solve(tol=1e-6, max_iter=200, verbose=True)
result.print_report()

state = NetworkState(
    fluid=water, node_P=result.node_P, node_h=result.node_h,
    node_mdot=result.node_mdot, params=result.params,
)

# Cycle corner states (1 = pump inlet, 2 = pump outlet / boiler inlet,
# 3 = boiler outlet / turbine inlet, 4 = turbine outlet / condenser inlet).
h1 = result.node_h["pump.in"]
h2 = result.node_h["pump.out"]
h3 = result.node_h["boiler.out"]
h4 = result.node_h["turbine.out"]

w_turbine = mdot * (h3 - h4)
w_pump = mdot * (h2 - h1)
q_in = mdot * (h3 - h2)
eta_cycle = (w_turbine - w_pump) / q_in

turb_metrics = turbine.report_metrics(state)

print("\nRankine cycle summary")
print(f"  turbine work    : {w_turbine / 1e6:8.3f} MW")
print(f"  pump work       : {w_pump / 1e3:8.1f} kW")
print(f"  boiler heat in  : {q_in / 1e6:8.3f} MW")
print(f"  net power       : {(w_turbine - w_pump) / 1e6:8.3f} MW")
print(f"  cycle efficiency: {eta_cycle * 100:7.2f} %")
print(f"  turbine exhaust quality x_out = {turb_metrics['x_out [-]']:.3f} (wet steam)")

print(
    "\nThis completes the Rankine cycle series. A real boiler that produces "
    "this steam has its own dynamics worth modeling on their own — see "
    "tutorials/06_boiler_drum_transient/ for the steam drum that would "
    "normally sit upstream of the boiler here."
)
