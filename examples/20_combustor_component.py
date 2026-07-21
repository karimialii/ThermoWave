"""SimpleCombustor: real fuel-mass-addition combustion, replacing the
"heater = Pipe with negative heat_loss" placeholder a plain heat-addition
cycle would otherwise use.

Single-shaft turboshaft: Source -> Compressor -> SimpleCombustor -> Turbine
-> Sink. Unlike the earlier shaft-coupling examples (which used a Controller
to find whatever compressor N hits a target exhaust temperature), this one
reflects a more typical real control split: shaft speed N is given directly
for both machines (e.g. held by a generator/grid frequency, so no free-N
unknown and no Shaft needed here — see 18_shaft_coupling.py for the case
where N itself is the solved unknown), and a Controller instead adjusts the
combustor's fuel flow (mdot_fuel left None, a free unknown) until a Sensor
on the turbine outlet reads the target exhaust temperature — i.e. fuel flow
is the actual temperature-control actuator, N is not. The heat addition
itself now really adds fuel mass flow (mdot_out = mdot_in + mdot_fuel) and
models a combustor pressure loss (PR < 1), rather than a massless enthalpy
bump.

A Cantera chemistry-based Combustor (equilibrium combustion products instead
of a fixed LHV) is available as a drop-in alternative when the optional
'cantera' extra is installed — see Combustor's docstring for what it does
differently and why that's a genuinely more accurate T_out, not just a
different API. See tutorials/02_gas_turbine_cycle/ for this same combustor
built into a full recuperated cycle with shaft coupling.

Run: .venv/bin/python examples/20_combustor_component.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.controller import Controller
from thermowave.components.sensor import Sensor
from thermowave.components.simple_combustor import SimpleCombustor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

shaft_N = 70000.0  # rev/min, held fixed for both machines (single rigid shaft)

src = Source(name="src", P=101325.0, T=288.15, mdot=0.8)
comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=gamma, N=shaft_N)
combustor = SimpleCombustor(name="cc", LHV=50.0e6, PR=0.96, efficiency=0.99, mdot_fuel=None)
turb = Turbine(name="turb", map_path="T100 Turb.tur", gamma=gamma, N=shaft_N)
snk = Sink(name="snk")

turb_outlet_sensor = Sensor(name="turb_outlet_sensor")
target_exhaust_T = 630 + 273.15
ctrl = Controller(
    name="ctrl",
    sensor=turb_outlet_sensor,
    quantity="T [K]",
    component=combustor,
    free_param="mdot_fuel",
    value=target_exhaust_T,
)

network = Network(fluid=air)
for component in (src, comp, combustor, turb, turb_outlet_sensor, ctrl, snk):
    network.add_component(component)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", combustor, "in")
network.connect(combustor, "out", turb, "in")
network.connect(turb, "out", turb_outlet_sensor, "tap")
network.connect(turb, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=400, damping=0.3, verbose=True)
result.print_report()

print(f"\nFuel flow: mdot_fuel = {result.params['cc.mdot_fuel'] * 1000.0:.2f} g/s")
