"""Sensor + Controller: a closed control loop, measurement decoupled from
actuation.

Setpoint (see 14_compressor_all_cases.py) ties a component's own
report_metrics() output straight to a target — e.g. a compressor's own
power or PR. Sensor + Controller generalizes that one step further: a Sensor
is a passive tap that can sit anywhere in the network reading (P, T, h,
mdot), and a Controller drives some *other* component's free parameter until
that independent reading hits a target — mirroring how a real plant's
measured point (a downstream thermocouple) and actuated point (a valve, a
shaft speed) are often different pieces of equipment.

This example drives a Compressor's shaft speed N until a Sensor sitting on
the compressor's own outlet reads a target temperature of 420 K.

Run: .venv/bin/python examples/17_sensor_controller.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.controller import Controller
from thermowave.components.sensor import Sensor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=gamma, N=None)
outlet_sensor = Sensor(name="outlet_sensor")
target_T = 420.0
ctrl = Controller(
    name="ctrl",
    sensor=outlet_sensor,
    quantity="T [K]",
    component=comp,
    free_param="N",
    value=target_T,
)
snk = Sink(name="snk")

network = Network(fluid=air)
for component in (src, comp, outlet_sensor, ctrl, snk):
    network.add_component(component)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", outlet_sensor, "tap")
network.connect(comp, "out", snk, "in")

result = network.solve(tol=1e-8, max_iter=200, damping=0.5, verbose=True)
result.print_report()

print(f"\nTarget outlet T = {target_T} K, controlled N = {result.params['comp.N']:.1f} rev/min")
