# Examples

Short, self-contained snippets — each one runnable as-is. For a full
end-to-end build with a thermal network and transient control, see the
[gas-turbine tutorial](../tutorials/building-a-gas-turbine-model.md).

## A single flow branch

The smallest possible network: a `Source` fixes inlet conditions, a `Pipe`
drops pressure by friction, a `Sink` terminates the branch.

```python
from thermowave.fluids.ideal_gas import IdealGasFluid
from thermowave.core.network import Network
from thermowave.components import Source, Pipe, Sink

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
net = Network(fluid=air)

src = Source(name="src", P=200_000.0, T=300.0, mdot=1.0)
pipe = Pipe(name="pipe", L=10.0, D=0.1, roughness=4.5e-5, mu=1.8e-5)
sink = Sink(name="sink")

for c in (src, pipe, sink):
    net.add_component(c)
net.connect(src, "out", pipe, "in")
net.connect(pipe, "out", sink, "in")

result = net.solve()
result.print_report()
```

## A throttling valve with a target downstream pressure

Give the valve's `opening` instead of a target, or flip it around: leave
`Sink(P=...)` fixed and let a `Setpoint` solve for the opening that hits it.
See [Control & instrumentation](../components/control/index.md) for the
`Setpoint`/`Controller` pattern this relies on.

```python
from thermowave.components import Source, Valve, Sink

src = Source(name="src", P=500_000.0, T=350.0, mdot=0.5)
valve = Valve(name="valve", D=0.05, K=8.0, opening=0.6)
sink = Sink(name="sink")

for c in (src, valve, sink):
    net.add_component(c)
net.connect(src, "out", valve, "in")
net.connect(valve, "out", sink, "in")

result = net.solve()
print(valve.report_metrics(result.state()))
```

## A recuperated Brayton cycle (analytic components)

The same shape as the [gas-turbine tutorial](../tutorials/building-a-gas-turbine-model.md),
using the fixed-PR/efficiency analytic components (`SimpleCompressor`,
`SimpleTurbine`, `SimpleCombustor`) instead of characteristic maps — a good
starting point before a real machine map is available.

```python
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid
from thermowave.components import (
    Source, SimpleCompressor, SimpleHeatExchanger, SimpleCombustor,
    SimpleTurbine, Sink, Setpoint,
)

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
net = Network(fluid=air)

src = Source(name="src", P=101_325.0, T=288.15, mdot=1.0)
comp = SimpleCompressor(name="comp", PR=4.0, eta_s=0.80)
recup = SimpleHeatExchanger(name="recup", effectiveness=0.80, PR_hot=0.98, PR_cold=0.98)
comb = SimpleCombustor(name="comb", LHV=50e6, mdot_fuel=None)
turb = SimpleTurbine(name="turb", PR=3.8, eta_s=0.85)
snk = Sink(name="snk", P=101_325.0)

for c in (src, comp, recup, comb, turb, snk):
    net.add_component(c)

net.connect(src, "out", comp, "in")
net.connect(comp, "out", recup, "cold_in")
net.connect(recup, "cold_out", comb, "in")
net.connect(comb, "out", turb, "in")
net.connect(turb, "out", recup, "hot_in")
net.connect(recup, "hot_out", snk, "in")

# Drive fuel flow to hit a target turbine inlet temperature.
target_tit = Setpoint(
    name="target_tit", component=comb, free_param="mdot_fuel",
    target_metric="T_out [K]", value=1150.0,
)
net.add_component(target_tit)

result = net.solve()
result.print_report()
```

## A PID-controlled setpoint change

Swap a steady `Controller` for a `PIDController` and step a setpoint —
see [Control & instrumentation](../components/control/index.md) for the full
contract (anti-windup, feedforward, `Schedule`-driven profiles).

```python
from thermowave.components import Sensor, PIDController

tit_sensor = Sensor(name="tit_sensor")
net.add_component(tit_sensor)
net.connect(comb, "out", tit_sensor, "tap")

pid = PIDController(
    name="pid_fuel", sensor=tit_sensor, quantity="T [K]",
    component=comb, free_param="mdot_fuel", setpoint=1150.0,
    Kp=1.0e-7, Ki=5.0e-8, Kd=0.0,
    output0=result.params["comb.mdot_fuel"],
)
net.replace_component(target_tit, pid)
history = net.solve_transient(duration=60.0, dt=0.5, initial=result)
history.plot((tit_sensor, "T [K]"))
```
