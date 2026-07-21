"""Tutorial 4, step 3: driving the power loop's setpoint through a profile.

Step 2 held pid_power's target fixed at 155 kW for the whole run. A real
plant's power target isn't fixed — it follows grid demand, changing over
time: hold, ramp, hold, step, hold. Rewriting the loop that changes it as a
hand-rolled Python for-loop around solve_transient() (manually reassigning
pid_power.setpoint at every step) would work, but Schedule does exactly
that: it owns a piecewise (t, value) profile and writes it into
pid_power.setpoint once per timestep, the same step() hook mechanism
PIDController itself uses to update its own output — so the whole run below
is still one solve_transient() call, not a loop you write by hand.

Profile: hold at 140 kW, ramp up to 155 kW, hold, ramp back down — small
enough excursions to stay inside the compressor-speed range this fixed-mdot
cycle behaves smoothly across (see step 2's note on why the power loop's
gain sign is what it is).

What you'll learn:
  - Schedule(target, attr, breakpoints) writes breakpoints[i][1] into
    target.attr at time breakpoints[i][0], holding flat before the first
    and after the last breakpoint, linearly interpolating between.
  - A Schedule is just another network component with a step() hook — same
    mechanism as PIDController and (from tutorial 3) a dynamic Shaft; the
    solver doesn't distinguish between them.
  - history.plot() can overlay the commanded setpoint against the actual
    response by reading schedule.value_at(t) directly, since the setpoint
    itself isn't a component metric.

Run: .venv/bin/python tutorials/04_closed_loop_control/03_power_profile_schedule.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.pid_controller import PIDController
from thermowave.components.schedule import Schedule
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
    setpoint=140000.0,   # overwritten every step by the Schedule below
    Kp=-0.02,
    Ki=-0.01,
    Kd=0.0,
    output0=58000.0,
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

# (t [s], target power [W]) breakpoints: hold at 140 kW, ramp to 155 kW,
# hold, ramp back down to 140 kW.
POWER_PROFILE_W = [
    (0.0, 140000.0), (10.0, 140000.0),
    (20.0, 155000.0),
    (30.0, 155000.0),
    (40.0, 140000.0),
    (50.0, 140000.0),
]
sched_power = Schedule(name="sched_power", target=pid_power, attr="setpoint", breakpoints=POWER_PROFILE_W)
network.add_component(sched_power)

history = network.solve_transient(duration=50.0, dt=1.0, tol=1e-8, max_iter=400, damping=0.3)

print(f"{'t [s]':>6}  {'P_set [kW]':>10}  {'P [kW]':>8}  {'N [rpm]':>8}  {'TOT [K]':>8}")
for i in range(0, len(history.times), 5):
    t, step = history.times[i], history.steps[i]
    state = NetworkState(
        fluid=step.fluid, node_P=step.node_P, node_h=step.node_h,
        node_mdot=step.node_mdot, params=step.params,
    )
    power_kw = gen.report_metrics(state)["power [W]"] / 1000.0
    tot = tot_sensor.report_metrics(state)["T [K]"]
    print(
        f"{t:6.1f}  {sched_power.value_at(t) / 1000.0:10.1f}  {power_kw:8.2f}  "
        f"{step.params['comp.N']:8.0f}  {tot:8.2f}"
    )

# Plotting requires the 'plot' extra: pip install thermowave[plot]. show=False
# + save_path so this stays non-interactive when run headless/in CI.
history.plot(
    (gen, "power [W]"), ylabel="power [W]", title="Electrical power vs. setpoint",
    show=False, save_path="power_profile_response.png",
)
print("\nSaved power_profile_response.png")

print(
    "\nThis closed-loop pattern — several PID loops plus a Schedule driving "
    "one of their setpoints through a profile, all inside one "
    "solve_transient() call — is the same structure a full production "
    "control-system model would use, just with a shorter, gentler profile."
)
