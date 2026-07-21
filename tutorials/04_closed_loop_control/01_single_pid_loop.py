"""Tutorial 4, step 1: one PID loop holding turbine outlet temperature constant.

Tutorial 2's cycle used a Setpoint — an exact, infinite-gain "controller"
that pins a free parameter to its target on every single solve, with no
time dimension at all. PIDController is the real time-domain counterpart:
inside Network.solve_transient(), it holds its actuated parameter fixed
during each algebraic solve, then updates it once per timestep from
Kp*error + Ki*integral + Kd*derivative — so it takes several steps to settle,
the way an actual control loop does, rather than jumping to the answer.

Here a PIDController reads a Sensor on the turbine outlet and actuates the
combustor's fuel flow (mdot_fuel) to hold turbine outlet temperature (TOT)
at a fixed 750 K target — more fuel means a hotter, more energetic
combustion product, so this is a direct, monotonic relationship, a natural
first loop to build (step 2 tackles a less direct one). It's started at a
fuel flow that gives a cooler-than-target TOT on purpose, so you can watch
it climb toward 750 K instead of already being there. Compressor and
turbine speed are both held fixed for now — tying them together on one
shaft, and controlling shaft speed itself, is step 2's job.

What you'll learn:
  - PIDController's sensor argument accepts any Sensor-shaped reading (see
    Sensor's own SENSOR_QUANTITIES) — here, "T [K]" off a turbine-outlet tap.
  - solve_transient() drives PIDController's step() hook once per timestep
    automatically — the same mechanism that integrates differential state
    (tutorial 3) also discovers and steps any step()-able component.
  - output_min/output_max keep the actuated fuel flow inside a physically
    sane range while the loop is still converging.

Run: .venv/bin/python tutorials/04_closed_loop_control/01_single_pid_loop.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.pid_controller import PIDController
from thermowave.components.sensor import Sensor
from thermowave.components.simple_combustor import SimpleCombustor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

shaft_N = 65000.0  # rev/min, held fixed for both machines (no shaft coupling yet)

src = Source(name="src", P=101325.0, T=288.15, mdot=0.8)
comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=gamma, N=shaft_N)
combustor = SimpleCombustor(name="cc", LHV=50.0e6, PR=0.96, efficiency=0.99, mdot_fuel=None)
turb = Turbine(name="turb", map_path="T100 Turb.tur", gamma=gamma, N=shaft_N)
tot_sensor = Sensor(name="tot_sensor")
snk = Sink(name="snk")

pid_fuel = PIDController(
    name="pid_fuel",
    sensor=tot_sensor,
    quantity="T [K]",
    component=combustor,
    free_param="mdot_fuel",
    setpoint=750.0,      # K target
    Kp=3.0e-6,
    Ki=1.0e-6,
    Kd=0.0,
    output0=0.008,       # kg/s, deliberately gives a cooler-than-target TOT
    output_min=0.004,
    output_max=0.020,
)

network = Network(fluid=air)
for component in (src, comp, combustor, turb, tot_sensor, pid_fuel, snk):
    network.add_component(component)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", combustor, "in")
network.connect(combustor, "out", turb, "in")
network.connect(turb, "out", tot_sensor, "tap")
network.connect(turb, "out", snk, "in")

history = network.solve_transient(duration=40.0, dt=1.0, tol=1e-8, max_iter=400, damping=0.4)

print(f"{'t [s]':>7}  {'mdot_fuel [kg/s]':>16}  {'TOT [K]':>8}")
for t, step in zip(history.times, history.steps):
    state = NetworkState(
        fluid=step.fluid, node_P=step.node_P, node_h=step.node_h,
        node_mdot=step.node_mdot, params=step.params,
    )
    tot = tot_sensor.report_metrics(state)["T [K]"]
    print(f"{t:7.1f}  {step.params['cc.mdot_fuel']:16.5f}  {tot:8.2f}")

print(
    "\nStep 2 adds a second, independent PID loop: shaft speed (now genuinely "
    "shared via Shaft) controlling electrical power, running at the same "
    "time as this fuel/TOT loop."
)
