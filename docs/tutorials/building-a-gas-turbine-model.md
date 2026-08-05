# Building a gas-turbine model

A recipe for a single-shaft recuperated microturbine (T100 class) in ThermoWave:
build it once, solve it steady, then run the transient from that steady state
with real controllers in the loop.

The shape is always the same seven steps:

1. create the network
2. add the flow-path components
3. add the thermal components
4. add the connections — flow, then heat
5. solve steady
6. swap the steady controllers for PID controllers
7. run the transient, seeded from the steady result

A runnable version of this against the repo's generic maps is
`tests/test_gas_turbine_workflow.py`, which follows these steps verbatim so the
guide can't drift from working code.

> **Maps and calibration constants.** Every number below is a placeholder. Real
> machine maps (`.cop`/`.tur`) and rig-fitted constants are yours to supply —
> nothing machine-specific ships in this repository.

---

## 1. Create the network

```python
from thermowave.core import Network
from thermowave.fluids import CanteraFluid

air = CanteraFluid(name="air", mechanism="gri30.yaml", composition="O2:0.21,N2:0.79")
network = Network(fluid=air)
```

Use a `CanteraFluid` if you want the combustor's reacted products to propagate
downstream (see *Composition-aware fluid propagation* in the README). An
`IdealGasFluid` is fine for a first cut and much faster.

## 2. Flow-path components

```python
from thermowave.components import (
    Combustor, Compressor, HeatExchanger, Sensor, Sink, Source, Turbine,
)

GAMMA = 1.4
P_AMB, T_AMB = 101325.0, 288.15

src   = Source(name="src", P=P_AMB, T=T_AMB, mdot=None, mdot_guess=0.7)
comp  = Compressor(name="comp", map_path="your-compressor.cop", gamma=GAMMA, N=None)
recup = HeatExchanger(name="recup", effectiveness=0.85, PR_hot=0.98, PR_cold=0.97)
comb  = Combustor(name="comb", PR=0.96, mdot_fuel=None, fuel="CH4")
turb  = Turbine(name="turb", map_path="your-turbine.tur", gamma=GAMMA, N=None)
tot   = Sensor(name="tot")          # turbine outlet temperature tap
snk   = Sink(name="snk", P=P_AMB)
```

**`None` is the "solve for this" sentinel.** `mdot=None`, `N=None`,
`mdot_fuel=None` each turn that quantity into a Newton unknown, and each then
needs exactly one equation to close it (step 5).

**`Shaft` and `PIDController` no longer read a component's free parameters at
construction time** — mechanical coupling is real port wiring now (`N` lives
on each turbomachine's own `"shaft"` port), so build order doesn't matter the
way it used to; what matters is wiring every mechanical/signal port
explicitly after `add_component()` (below).

## 3. Shaft and load

```python
from thermowave.components import Shaft, ShaftLoad

ETA_GEN, ETA_MECH = 0.96, 0.98

load  = ShaftLoad(name="load", power=75_000.0, efficiency=ETA_GEN * ETA_MECH)
shaft = Shaft(
    name="shaft", members=[comp, turb, load],
    efficiency=1.0, inertia=0.05, dynamic=True, N0=62_000.0,
)
```

`signs` is left at its default (`None`) here — it's inferred from each
member's own `shaft_sign()`: turbine `+1` (delivers), compressor `-1`
(absorbs), load `-1` (draws), the same answer you'd get writing
`signs=[-1.0, 1.0, -1.0]` by hand. Pass `signs` explicitly only to override
one of these for a genuine edge case (a motor running in regenerative mode,
say). `gear_ratios`, unlike `signs`, has one entry per *speed-tied* member
only (`load` has no `"shaft"` mechanical port, so it's torque-only and
doesn't count there).

`Shaft` exposes one mechanical port (`"m0"`, `"m1"`, ...) per speed-tied
member and one signal port (`"p0"`, `"p1"`, ...) per member, in `members`
order. Connecting a member's mechanical port pulls its power signal in
automatically — a real shaft doesn't have a separate wire for power either —
so `comp`/`turb` only need one `connect()` call each; `load` has no
`"shaft"` mechanical port to piggyback on, so it still needs its own
explicit signal connection:

```python
network.connect(shaft, "m0", comp, "shaft", kind="mechanical")
network.connect(shaft, "m1", turb, "shaft", kind="mechanical")
network.connect(shaft, "p2", load, "power", kind="signal")
```

(Or skip writing these by hand entirely: `shaft.autowire(network)`, called
once after `add_component(shaft)`, makes all three calls for you.)

Use `ShaftLoad`, not `SimpleGenerator`, on a dynamic shaft. `SimpleGenerator` is
a passive reader — it reports power but exerts no torque, so it never closes the
power loop. `ShaftLoad` is a real boundary condition: raise its `power` and the
equilibrium speed moves.

See the mechanical-loss trap below for why `efficiency` is on the load and
`Shaft(efficiency=1.0)`.

## 4a. Thermal components

```python
from thermowave.components import Conduction, Convection, ThermalMass

comp_casing = ThermalMass(name="comp_casing", thermal_capacitance=2000.0, T0=T_AMB)
turb_casing = ThermalMass(name="turb_casing", thermal_capacitance=2000.0, T0=T_AMB)
shaft_mass  = ThermalMass(name="shaft_mass",  thermal_capacitance=500.0,  T0=T_AMB)
```

`thermal_capacitance` is lumped `m·cp` [J/K]. In a steady solve each mass floats
to whatever temperature makes its net heat zero; in a transient it integrates
forward from `T0`.

## 4b. Connections

Flow first:

```python
for component in (src, comp, recup, comb, turb, tot, load, shaft, snk):
    network.add_component(component)

network.connect(src,   "out",      comp,  "in")
network.connect(comp,  "out",      recup, "cold_in")
network.connect(recup, "cold_out", comb,  "in")
network.connect(comb,  "out",      turb,  "in")
network.connect(turb,  "out",      tot,   "tap")
network.connect(turb,  "out",      recup, "hot_in")
network.connect(recup, "hot_out",  snk,   "in")
```

`turb.out` is connected twice on purpose — the sensor tap and the recuperator
hot side are the same physical node, and `connect()` merges them.

Then heat:

```python
network.add_heat_path(Convection(name="conv_comp", a=(comp, "out"), b=comp_casing, h=50.0, A=0.3))
network.add_heat_path(Convection(name="conv_turb", a=(turb, "out"), b=turb_casing, h=50.0, A=0.3))
network.add_heat_path(Conduction(name="cond_ts", a=turb_casing, b=shaft_mass, k=15.0, A=0.01, L=0.4))
network.add_heat_path(Conduction(name="cond_sc", a=shaft_mass, b=comp_casing, k=15.0, A=0.01, L=0.4))
network.add_heat_path(Convection(name="conv_ta", a=turb_casing, b=T_AMB, h=10.0, A=1.0))
network.add_heat_path(Convection(name="conv_ca", a=comp_casing, b=T_AMB, h=10.0, A=1.0))
```

An endpoint is a `ThermalMass`, a `(component, port)` pair for a fluid node, or a
bare float for ambient. `Q` is positive when `a` is hotter than `b`.

`add_heat_path()` registers the path and any `ThermalMass` endpoint, and derives
each endpoint's sign from whether it is `a` or `b`. Those signs are mirror
images — a mass *gains* as `b`, a flow component *loses* as `a` — which is
exactly the bookkeeping that used to be hand-written per path, unvalidated, and
silently wrong if flipped.

Before solving, ask:

```python
assert network.check_wiring() == []
```

It names any free parameter nothing closes, any closed by two components at
once, and any controller left pointing at a removed component.

## 5. Solve steady

A dynamic shaft closes its own speed (`state_derivative() == 0` at steady
state), so the remaining unknown is fuel:

```python
from thermowave.components import Controller

ctrl_tot = Controller(
    name="ctrl_tot", sensor=tot, quantity="T [K]",
    component=comb, free_param="mdot_fuel", value=918.15,   # 645 °C
)
network.add_component(ctrl_tot)

steady = network.solve(tol=1e-7, max_iter=800, damping=0.3)
```

A `Controller` is an ideal, infinite-gain controller: it hits the target exactly
on every solve, with no dynamics. That's what you want for a steady operating
point.

If a cold start won't converge, solve an easier version first and pass its
result as `warm_start=` — continuation, the standard fix for a tightly-coupled
turbomachinery network.

## 6. Swap in the PID controllers

```python
from thermowave.components import PIDController

pid_fuel = PIDController(
    name="pid_fuel", sensor=tot, quantity="T [K]",
    component=comb, free_param="mdot_fuel", setpoint=918.15,
    Kp=6.0e-6, Ki=3.0e-6, Kd=0.0,
    output0=steady.params["comb.mdot_fuel"],     # seed from the steady point
    output_min=0.0012, output_max=0.0080,
)
network.replace_component(ctrl_tot, pid_fuel)
```

`replace_component()` matters here. `Controller` and `PIDController` each
contribute one residual pinning `comb.mdot_fuel`; adding the PID without
removing the `Controller` leaves two equations for one unknown, and the solve
fails as non-square.

Seeding `output0` from the steady result is what makes step 7 start from rest
instead of lurching: `output0` is both the initial output *and* the permanent
bias the PID trims around.

## 7. Run the transient

```python
history = network.solve_transient(
    duration=200.0, dt=1.0, initial=steady,
    tol=1e-7, max_iter=3000, damping=0.15,
)
```

`initial=steady` supplies every differential state's `t=0` value — the shaft
speed and each casing temperature. Because this is the *same* network object the
steady solve ran on, those states line up by name automatically.

For a load profile, drive the load with a `Schedule`:

```python
from thermowave.components import Schedule

network.add_component(Schedule(
    name="load_profile", target=load, attr="power",
    breakpoints=[(0.0, 50e3), (60.0, 50e3), (90.0, 75e3), (200.0, 75e3)],
))
```

Anything with a `step(state, dt)` method is discovered and stepped once per
accepted timestep. Add `adaptive=True` if the run has both a fast initial
transient and a slow settle.

---

## Four traps worth knowing before you debug them

### Mechanical loss belongs on the load, not the shaft

`Shaft.efficiency` multiplies the **net** shaft power. On a dynamic shaft the
steady equilibrium is exactly net == 0, so the factor cancels and the loss
silently disappears. Measured on a T100 at 75 kW / 62,321 rpm: `ETA_MECH=0.98`
on the `Shaft` asked for 0.00561 kg/s of fuel where the correct build wants
0.00567 — the whole 2% loss missing.

Put it on the load: `ShaftLoad(efficiency=ETA_GEN * ETA_MECH)` with
`Shaft(efficiency=1.0)`. The shaft-side draw becomes `P_elec / (ETA_GEN ·
ETA_MECH)`, so the turbine surplus genuinely carries it. ThermoWave warns if you
set `Shaft(dynamic=True, efficiency != 1.0)`.

Also: don't read electrical power from `SimpleGenerator` in dynamic mode. It
reads the shaft's net power, which is ~0 W at equilibrium. Electrical power is
`ShaftLoad.power` — it's the commanded demand, so it's exact — and the shaft's
net power becomes a torque-imbalance diagnostic that should approach zero.

### More fuel *lowers* turbine outlet temperature at fixed load

Measured at the 75 kW point: −1% fuel → 698 °C, nominal → 645 °C, +1% fuel →
573 °C. At *fixed load* the extra energy goes into shaft speed rather than
output, and the resulting extra air flow dilutes the temperature rise faster
than the extra fuel raises it.

So TOT is not a monotone "more fuel is hotter" knob, and a fuel→TOT control
pairing has the sign you would not guess. Check the sign of `Kp` against the
plant before blaming the gains.

### There is a no-equilibrium cliff just below the operating fuel flow

At −2% fuel there is **no steady state at all** — the shaft cannot sustain a
constant 75 kW draw on that little fuel. A Newton solve pushed there walks to
negative temperature and raises `CanteraError` from the combustor's equilibrium
call, **not** `ConvergenceError`. Solve wrappers must catch `Exception`.

This is what kills a naive PID: with no feedforward, the integral term has to
wander the whole fuel range to find the operating point, and wandering ~2% low
falls off the cliff. It is not a matter of gains.

The structural fix is to bias the governor from the **load**:

```python
pid_fuel = PIDController(
    ...,
    feedforward=lambda state: fuel_schedule(load.power),
)
```

`feedforward` replaces the constant bias with one computed fresh each step, so
the PI only ever trims a small correction around a known-good operating point.
Bias from the *load*, never from a rate-limited commanded speed: measured
failure of the latter — during a load ramp it fed fuel for ~52 kW while the load
was already 55.8 kW, so the shaft decelerated instead of accelerating, TOT ran
641 → 726 °C, and the solve diverged. The fuel a machine needs is set by the
energy balance, i.e. by the load it is carrying.

### Sweeps must run ascending

A descending cold start fails at 100 kW with a singular Jacobian. Walk power
upward and thread each result into the next as `warm_start=`.

---

## Reading results

```python
state = steady.state()
print(comp.report_metrics(state)["PR [-]"])
print(shaft.report_metrics(state)["power [W]"])   # ~0 at dynamic equilibrium
print(load.power)                                 # the real electrical output

steady.print_report()
history.plot((tot, "T [K]"), (shaft, "N [rev/min]"))
```

`SolveResult.state()` rebuilds the `NetworkState` any `report_metrics()` call
needs — a converged result stores node values, but anything a component
*computes* (pressure ratio, heat flux, net power) has to be evaluated against a
state.
