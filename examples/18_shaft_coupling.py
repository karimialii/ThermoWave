"""Shaft: mechanical connection (with inertia and efficiency) between two
or more components sharing one physical shaft.

A single-shaft turboshaft: the Compressor and Turbine no longer get
independent N values (fixed or via their own separate Setpoint) — instead
both are left fully free (N=None), a Shaft ties comp.N == turb.N
(gear_ratios=[1.0], direct coupling), and a single Controller pins down
that now-shared speed by targeting a Sensor's reading of the turbine outlet
temperature. This mirrors how a real single-shaft engine is actually
constrained: one physical shaft speed, set by whatever the control system
is actually regulating (here, exhaust temp), not two independently guessed
numbers.

signs=[-1.0, 1.0] tells the Shaft the compressor draws power from it while
the turbine delivers power to it, so its own report_metrics()["power [W]"]
reads as the net mechanical power actually available on the shaft (turbine
output minus compressor draw), scaled down by efficiency=0.98 for
bearing/windage losses. inertia is set but — this is Shaft's default
dynamic=False (steady-state, control-loop) mode, so it has no effect here.
See tutorials/03_transient_rotor_dynamics/ for dynamic=True, where Shaft
instead owns its own speed as differential state driven by inertia and net
torque, with no external Controller/Setpoint needed for the shaft at all —
and tutorials/02_gas_turbine_cycle/ for this same shaft-coupling pattern
built up step by step into a full recuperated cycle.

Run: .venv/bin/python examples/18_shaft_coupling.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.controller import Controller
from thermowave.components.pipe import Pipe
from thermowave.components.sensor import Sensor
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
    name="shaft",
    components=[comp, turb],
    signs=[-1.0, 1.0],
    efficiency=0.98,
    inertia=0.05,
)
snk = Sink(name="snk")

turb_outlet_sensor = Sensor(name="turb_outlet_sensor")
target_exhaust_T = 900.0
ctrl = Controller(
    name="ctrl",
    sensor=turb_outlet_sensor,
    quantity="T [K]",
    component=comp,
    free_param="N",
    value=target_exhaust_T,
)

network = Network(fluid=air)
for component in (src, comp, heater, turb, shaft, turb_outlet_sensor, ctrl, snk):
    network.add_component(component)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", heater, "in")
network.connect(heater, "out", turb, "in")
network.connect(turb, "out", turb_outlet_sensor, "tap")
network.connect(turb, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=400, damping=0.3, verbose=True)
result.print_report()

print(
    f"\nShared shaft speed: comp.N = {result.params['comp.N']:.1f} rev/min, "
    f"turb.N = {result.params['turb.N']:.1f} rev/min"
)
print(f"Target exhaust T = {target_exhaust_T} K")
