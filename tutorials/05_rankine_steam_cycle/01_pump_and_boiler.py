"""Tutorial 5, step 1: pressurizing and boiling water.

Every component tutorial so far has used a single-phase ideal gas. This
tutorial switches to CoolPropFluid (real water/steam properties, including
the two-phase dome) and the phase-change components: a Pump raises
feedwater pressure, and a SimpleEvaporator boils it into superheated steam.

Pump is entropy-based (works for liquids and wet steam, unlike the gas
turbines' gamma relation): it computes the ideal (isentropic) pressure rise
via the fluid's own entropy_ph()/enthalpy_ps() methods, then scales the
actual enthalpy rise up by 1/eta for inefficiency. SimpleEvaporator adds
heat until the outlet hits a target condition — here, saturated vapor plus
150 K of superheat — rather than taking a duty as a direct input, since cp
is effectively infinite during a constant-pressure phase change (see its own
docstring for why an effectiveness/UA calculation can't represent boiling).

What you'll learn:
  - CoolPropFluid(name="Water") gives real substance properties; a Source's
    (P, T) can pin any single-phase state but not an exactly-saturated or
    two-phase one — feedwater here is specified slightly subcooled.
  - Pump(P_out=..., eta=...) raises pressure with a small enthalpy rise.
  - SimpleEvaporator(superheat=...) boils and superheats in one component,
    reporting the resulting quality/temperature rather than taking them as
    inputs.

Requires the 'coolprop' extra: pip install thermowave[coolprop]

Run: .venv/bin/python tutorials/05_rankine_steam_cycle/01_pump_and_boiler.py
"""

from thermowave.components.pump import Pump
from thermowave.components.simple_evaporator import SimpleEvaporator
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.fluids.real_fluid import CoolPropFluid

water = CoolPropFluid(name="Water")

P_cond = 1.0e4    # 0.1 bar
P_boiler = 4.0e6  # 40 bar

feed = Source(name="feed", P=P_cond, T=315.0, mdot=10.0)  # slightly subcooled liquid
pump = Pump(name="pump", P_out=P_boiler, eta=0.75)
boiler = SimpleEvaporator(name="boiler", superheat=150.0)
snk = Sink(name="snk")

network = Network(fluid=water)
for component in (feed, pump, boiler, snk):
    network.add_component(component)

network.connect(feed, "out", pump, "in")
network.connect(pump, "out", boiler, "in")
network.connect(boiler, "out", snk, "in")

result = network.solve(tol=1e-6, max_iter=200, verbose=True)
result.print_report()

print(
    "\nSteam now leaves the boiler at 40 bar, superheated — but nothing "
    "extracts work from it yet. Step 2 adds a SteamTurbine to expand it."
)
