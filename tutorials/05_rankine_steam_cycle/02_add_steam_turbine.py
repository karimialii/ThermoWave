"""Tutorial 5, step 2: expanding steam through a turbine.

Step 1's superheated steam now expands through a SteamTurbine back down to
condenser pressure, extracting shaft work. SteamTurbine is the wet-steam-
correct counterpart to the gas turbines' SimpleTurbine: it uses the same
entropy-based isentropic path Pump does (entropy_ph()/enthalpy_ps()) rather
than an ideal-gas temperature relation, which would be wrong once the
expansion crosses into the two-phase dome — exactly what happens here.

report_metrics()["x_out [-]"] is the turbine exhaust quality: how much of
the exit flow is still vapor (1.0) vs. condensed liquid (0.0). Real steam
turbines care about this a great deal — too much liquid in the last stages
erodes the blades — so it's reported directly rather than requiring you to
compute it from (P, h) yourself.

What you'll learn:
  - SteamTurbine(P_out=..., eta_s=...) mirrors SimpleTurbine's constructor
    shape (a pressure target and an isentropic efficiency), just built on
    entropy instead of gamma.
  - report_metrics()["x_out [-]"] tells you directly whether an expansion
    ended up in the wet region — no separate quality calculation needed.

Requires the 'coolprop' extra: pip install thermowave[coolprop]

Run: .venv/bin/python tutorials/05_rankine_steam_cycle/02_add_steam_turbine.py
"""

from thermowave.components.pump import Pump
from thermowave.components.simple_evaporator import SimpleEvaporator
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.steam_turbine import SteamTurbine
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.real_fluid import CoolPropFluid

water = CoolPropFluid(name="Water")

P_cond = 1.0e4
P_boiler = 4.0e6

feed = Source(name="feed", P=P_cond, T=315.0, mdot=10.0)
pump = Pump(name="pump", P_out=P_boiler, eta=0.75)
boiler = SimpleEvaporator(name="boiler", superheat=150.0)
turbine = SteamTurbine(name="turbine", P_out=P_cond, eta_s=0.85)
snk = Sink(name="snk")

network = Network(fluid=water)
for component in (feed, pump, boiler, turbine, snk):
    network.add_component(component)

network.connect(feed, "out", pump, "in")
network.connect(pump, "out", boiler, "in")
network.connect(boiler, "out", turbine, "in")
network.connect(turbine, "out", snk, "in")

result = network.solve(tol=1e-6, max_iter=200, verbose=True)
result.print_report()

state = NetworkState(
    fluid=water, node_P=result.node_P, node_h=result.node_h,
    node_mdot=result.node_mdot, params=result.params,
)
metrics = turbine.report_metrics(state)
print(f"\nTurbine power = {metrics['power [W]'] / 1e6:.3f} MW")
print(f"Turbine exhaust quality x_out = {metrics['x_out [-]']:.3f} (wet steam)")

print(
    "\nThe turbine exhaust leaves as a Sink here, discarding the low-grade "
    "heat still in that wet steam. Step 3 condenses it back to liquid "
    "instead, completing the cycle and closing the loop back to feedwater "
    "conditions — the point where a real cycle efficiency can be computed."
)
