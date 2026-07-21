"""Tutorial 2, step 3: adding a recuperator and a real fuel supply line.

Two upgrades over step 2's shaft-coupled cycle, both common in small
recuperated gas turbines like the one whose maps these examples use:

1. A recuperator (SimpleHeatExchanger) preheats the compressor discharge
   using heat recovered from the turbine exhaust before it's wasted —
   "cold_in"/"cold_out" is the compressor-discharge side being heated,
   "hot_in"/"hot_out" is the turbine-exhaust side giving up heat. This is
   the same SimpleHeatExchanger any two streams could use; nothing about it
   is turbine-specific.

2. Fuel now arrives through the combustor's own second inlet
   (SimpleCombustor(use_fuel_port=True)) fed by a real Source -> Pipe fuel
   line, instead of step 1/2's plain mdot_fuel number. That fuel branch has
   its own pressure, temperature, and line pressure drop — so fuel's own
   sensible enthalpy enters the energy balance too, not just its heating
   value, and the combustor's fuel_in port is just another node like any
   other component's inlet.

A Sensor reads the turbine outlet temperature (TOT) — a real engine's most
watched number, since it limits how hard the turbine can be run without
damaging the hot-section blades.

What you'll learn:
  - SimpleHeatExchanger's four named ports (hot_in/hot_out/cold_in/cold_out)
    connect two independent streams through one component.
  - A component can expose more than the default "in"/"out" ports —
    SimpleCombustor(use_fuel_port=True) adds "fuel_in".
  - Sensor is a zero-residual "tap": connecting it to a node just gives you
    a read of that node's state without perturbing the solve.

Run: .venv/bin/python tutorials/02_gas_turbine_cycle/03_recuperated_cycle.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.pipe import Pipe
from thermowave.components.sensor import Sensor
from thermowave.components.setpoint import Setpoint
from thermowave.components.shaft import Shaft
from thermowave.components.simple_combustor import SimpleCombustor
from thermowave.components.simple_heat_exchanger import SimpleHeatExchanger
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

target_N = 65000.0  # rev/min

src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=gamma, N=None)
recup = SimpleHeatExchanger(name="recup", effectiveness=0.8, PR_hot=0.98, PR_cold=0.97)

fuel_src = Source(name="fuel_src", P=350000.0, T=300.0, mdot=0.0087)
fuel_pipe = Pipe(name="fuel_pipe", L=1.0, D=0.01, f=0.02, n_elem=1)
combustor = SimpleCombustor(name="cc", LHV=50.0e6, PR=0.97, efficiency=0.99, use_fuel_port=True)

turb = Turbine(name="turb", map_path="T100 Turb.tur", gamma=gamma, N=None)
shaft = Shaft(name="shaft", components=[comp, turb], signs=[-1.0, 1.0], efficiency=0.98)
setpoint = Setpoint(
    name="sp_N", component=comp, free_param="N", target_metric="N [rev/min]", value=target_N
)
tot_sensor = Sensor(name="tot_sensor")
snk = Sink(name="snk")

network = Network(fluid=air)
for component in (
    src, comp, recup, fuel_src, fuel_pipe, combustor, turb, shaft, setpoint, tot_sensor, snk,
):
    network.add_component(component)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", recup, "cold_in")
network.connect(recup, "cold_out", combustor, "in")
network.connect(fuel_src, "out", fuel_pipe, "in")
network.connect(fuel_pipe, "out", combustor, "fuel_in")
network.connect(combustor, "out", turb, "in")
network.connect(turb, "out", tot_sensor, "tap")
network.connect(turb, "out", recup, "hot_in")
network.connect(recup, "hot_out", snk, "in")

result = network.solve(tol=1e-8, max_iter=100, damping=1.0, verbose=True)
result.print_report()

state = NetworkState(
    fluid=air, node_P=result.node_P, node_h=result.node_h,
    node_mdot=result.node_mdot, params=result.params,
)
print(f"\nTOT (turbine outlet temperature) = {tot_sensor.report_metrics(state)['T [K]']:.2f} K")

print(
    "\nThis is the steady-state cycle this whole tutorial series builds on:\n"
    "  - tutorials/03_transient_rotor_dynamics/ makes the shaft dynamic and\n"
    "    watches it spool up over time instead of assuming equilibrium.\n"
    "  - tutorials/04_closed_loop_control/ replaces this cycle's fixed fuel\n"
    "    flow and shaft-speed Setpoint with real PID control loops."
)
