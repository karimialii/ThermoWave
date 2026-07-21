"""PIDController: finite-response closed-loop control, vs. Controller's
instantaneous/ideal one.

Controller (see 17_sensor_controller.py) pins a Sensor reading to a target
exactly on every solve — a perfect, infinite-gain controller with no
dynamics. PIDController is the time-domain counterpart: inside
Network.solve_transient(), it holds a component's free parameter fixed
during each algebraic solve, then updates that value once per timestep from
Kp*error + Ki*integral + Kd*derivative — so it takes several steps to settle,
same as a real control loop.

Here a Compressor's shaft speed N (left free) is driven by a PID loop reading
a Sensor on the compressor outlet, targeting an outlet temperature of 420 K,
starting from an initial guess (N=60000 rev/min) that's actually too hot.
Watch the outlet temperature settle toward the 420 K target over the run
instead of jumping there in one step.

See tutorials/04_closed_loop_control/ for this same PID pattern applied to
a full gas turbine cycle, with multiple loops and a time-varying setpoint.

Run: .venv/bin/python examples/21_pid_controller.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.pid_controller import PIDController
from thermowave.components.sensor import Sensor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=gamma, N=None)
sensor = Sensor(name="outlet_sensor")
snk = Sink(name="snk")

target_T = 420.0
pid = PIDController(
    name="pid",
    sensor=sensor,
    quantity="T [K]",
    component=comp,
    free_param="N",
    setpoint=target_T,
    Kp=60.0,
    Ki=50.0,
    Kd=0.0,
    output0=60000.0,       # rev/min, deliberately too hot (T ~ 445 K there)
    output_min=10000.0,
    output_max=100000.0,
)

network = Network(fluid=air)
for component in (src, comp, sensor, pid, snk):
    network.add_component(component)

network.connect(src, "out", comp, "in")
network.connect(comp, "out", sensor, "tap")
network.connect(comp, "out", snk, "in")

history = network.solve_transient(
    duration=8.0, dt=0.1, tol=1e-6, max_iter=300, damping=0.5,
)

print(f"{'t [s]':>8}  {'N [rev/min]':>12}  {'T_out [K]':>10}")
for t, step in zip(history.times, history.steps):
    N = step.params["comp.N"]
    T_out = air.temperature_ph(step.node_P["comp.out"], step.node_h["comp.out"])
    print(f"{t:8.2f}  {N:12.1f}  {T_out:10.2f}")

print(f"\nTarget T_out = {target_T} K")

# Plotting requires the 'plot' extra: pip install thermowave[plot]
history.plot(
    (sensor, "T [K]"), ylabel="T [K]", title="PID-controlled compressor outlet temperature",
    show=False, save_path="pid_temperature.png",
)
print("Saved pid_temperature.png")
