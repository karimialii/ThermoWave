"""Generator + SimpleGenerator: electrical power from a turbine's shaft.

Extends 18_shaft_coupling.py's single-shaft turboshaft (Compressor + Turbine
locked to one shaft speed, pinned by a Controller reading turbine outlet
temperature) with two alternative generator models hung off the same
Turbine, both purely passive readers (no ports, no residuals — they don't
feed back into the thermodynamic solve):

- SimpleGenerator: takes the Turbine's own reported shaft power and scales
  it by a fixed mechanical-to-electrical efficiency.
- Generator: ignores the Turbine's reported power entirely and instead
  looks up torque at the Turbine's shaft speed from its own speed-vs-torque
  map ("Generator Torque.gen"), computing power_mech = torque * omega itself
  — useful when the generator's own rating curve is what's actually
  limiting/driving electrical output, not the turbine's thermodynamic power.

Run: .venv/bin/python examples/19_generator_component.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.controller import Controller
from thermowave.components.generator import Generator
from thermowave.components.pipe import Pipe
from thermowave.components.sensor import Sensor
from thermowave.components.shaft import Shaft
from thermowave.components.simple_generator import SimpleGenerator
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
shaft = Shaft(name="shaft", components=[comp, turb], signs=[-1.0, 1.0], efficiency=0.98)
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

simple_gen = SimpleGenerator(name="simple_gen", component=turb, efficiency=0.95)
map_gen = Generator(name="map_gen", component=turb, map_path="Generator Torque.gen", efficiency=0.95)

network = Network(fluid=air)
for component in (
    src,
    comp,
    heater,
    turb,
    shaft,
    turb_outlet_sensor,
    ctrl,
    simple_gen,
    map_gen,
    snk,
):
    network.add_component(component)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", heater, "in")
network.connect(heater, "out", turb, "in")
network.connect(turb, "out", turb_outlet_sensor, "tap")
network.connect(turb, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=400, damping=0.3, verbose=True)
result.print_report()

print(
    f"\nTurbine shaft power = {result.params.get('turb.N', float('nan')):.1f} rev/min shaft speed"
)
