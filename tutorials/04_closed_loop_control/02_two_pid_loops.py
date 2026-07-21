"""Tutorial 4, step 2: two independent PID loops running at once.

Step 1's fuel/TOT loop returns unchanged. This adds a second, completely
independent loop: compressor and turbine are now genuinely tied together on
one Shaft (as in tutorial 2), and a second PIDController reads
SimpleGenerator's electrical power output and actuates the shared shaft
speed N to hold power at a target.

Two PID loops wired into the same network, with no interaction between them
beyond the physics itself, is exactly how real closed-loop gas turbine
control works: each loop reads its own sensor and drives its own actuator,
and solve_transient() steps every step()-able component once per timestep
regardless of how many there are.

One sign worth noticing: this loop's Kp and Ki are *negative*. In this
particular network, mass flow is held fixed at the Source (mdot=0.8 kg/s,
not left free) — so raising N raises the compressor's pressure ratio without
also pulling more air through, which in this operating region actually
*reduces* net shaft power rather than increasing it. The controller has to
push N in the opposite direction from what "more speed, more power"
intuition would suggest, precisely because of that fixed-mass-flow choice
— a concrete reminder that a PID loop's sign always follows the actual
local sensitivity of output to input, not a general rule of thumb.

What you'll learn:
  - Multiple PIDControllers coexist in one network with zero special wiring
    — each is just another component with its own step() hook.
  - A control loop's gain sign is a property of the actual system it's
    controlling, not a fixed convention — always check the real sensitivity
    before assuming the sign.

Run: .venv/bin/python tutorials/04_closed_loop_control/02_two_pid_loops.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.pid_controller import PIDController
from thermowave.components.sensor import Sensor
from thermowave.components.shaft import Shaft
from thermowave.components.simple_combustor import SimpleCombustor
from thermowave.components.simple_generator import SimpleGenerator
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

src = Source(name="src", P=101325.0, T=288.15, mdot=0.8)
comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=gamma, N=None)
combustor = SimpleCombustor(name="cc", LHV=50.0e6, PR=0.96, efficiency=0.99, mdot_fuel=None)
turb = Turbine(name="turb", map_path="T100 Turb.tur", gamma=gamma, N=None)
shaft = Shaft(name="shaft", components=[comp, turb], signs=[-1.0, 1.0], efficiency=0.98)
gen = SimpleGenerator(name="gen", component=shaft, efficiency=0.95)
tot_sensor = Sensor(name="tot_sensor")
snk = Sink(name="snk")

pid_fuel = PIDController(
    name="pid_fuel",
    sensor=tot_sensor,
    quantity="T [K]",
    component=combustor,
    free_param="mdot_fuel",
    setpoint=750.0,
    Kp=3.0e-6,
    Ki=1.0e-6,
    Kd=0.0,
    output0=0.010,
    output_min=0.004,
    output_max=0.020,
)
pid_power = PIDController(
    name="pid_power",
    sensor=gen,
    quantity="power [W]",
    component=comp,
    free_param="N",
    setpoint=155000.0,   # 155 kW target
    Kp=-0.02,
    Ki=-0.01,
    Kd=0.0,
    output0=58000.0,     # rev/min
    output_min=50000.0,
    output_max=65000.0,
)

network = Network(fluid=air)
for component in (
    src, comp, combustor, turb, shaft, gen, tot_sensor, pid_fuel, pid_power, snk,
):
    network.add_component(component)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", combustor, "in")
network.connect(combustor, "out", turb, "in")
network.connect(turb, "out", tot_sensor, "tap")
network.connect(turb, "out", snk, "in")

history = network.solve_transient(duration=20.0, dt=1.0, tol=1e-8, max_iter=400, damping=0.3)

print(f"{'t [s]':>6}  {'N [rpm]':>8}  {'mdot_fuel':>10}  {'P [kW]':>8}  {'TOT [K]':>8}")
for t, step in zip(history.times, history.steps):
    state = NetworkState(
        fluid=step.fluid, node_P=step.node_P, node_h=step.node_h,
        node_mdot=step.node_mdot, params=step.params,
    )
    power_kw = gen.report_metrics(state)["power [W]"] / 1000.0
    tot = tot_sensor.report_metrics(state)["T [K]"]
    print(
        f"{t:6.1f}  {step.params['comp.N']:8.0f}  {step.params['cc.mdot_fuel']:10.5f}  "
        f"{power_kw:8.2f}  {tot:8.2f}"
    )

print(
    "\nBoth targets (155 kW, 750 K) are climbed toward independently. Step 3 "
    "drives pid_power's own setpoint through a time-varying profile with a "
    "Schedule component, instead of holding one fixed target."
)
