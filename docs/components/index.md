# Components

Every component in ThermoWave implements the same interface, `BaseComponent`
(`thermowave.components.base_component`). Understanding that interface makes
every component page below read the same way: what its **ports** are, what
**equations** it contributes to the Newton system, and what (if anything) it
leaves **free** for something else in the network to pin down.

```{toctree}
:maxdepth: 2

flow-elements/index
turbomachinery/index
mechanical/index
combustion/index
heat-exchangers/index
thermal-network/index
control/index
```

| Category | Components |
|---|---|
| [Flow elements](flow-elements/index.md) | `Source`, `Sink`, `Pipe`, `Valve`, `CheckValve`, `Nozzle`, `Junction` |
| [Turbomachinery](turbomachinery/index.md) | `SimpleCompressor`, `Compressor`, `SimpleTurbine`, `Turbine`, `SteamTurbine`, `Pump` |
| [Mechanical & electrical](mechanical/index.md) | `Shaft`, `ShaftLoad`, `ElectricMotor`, `Generator`, `SimpleGenerator` |
| [Combustion](combustion/index.md) | `SimpleCombustor`, `Combustor` |
| [Heat exchangers & phase change](heat-exchangers/index.md) | `SimpleHeatExchanger`, `MultiPassHeatExchanger`, `Condenser`, `SimpleCondenser`, `Evaporator`, `SimpleEvaporator`, `Drum`, `Tank` |
| [Thermal network](thermal-network/index.md) | `ThermalMass`, `Convection`, `Conduction`, `Radiation` |
| [Control & instrumentation](control/index.md) | `Sensor`, `Controller`, `PIDController`, `Setpoint`, `Schedule` |

Each category page lists its components with a one-line description; each
component then has its own page with its diagram, ports, and equations.

## The modeling framework

### Nodes, ports, and connections

A `Network` holds the component list and a union-find over **node ids**.
Every component names its own local ports in `ports() -> dict[str, str]`
(conventionally `f"{self.name}.{port}"`, e.g. `"comp.in"`). `Network.connect(a,
"out", b, "in")` merges two components' local port ids into one shared node
— after that, both components read and write the same `(P, h)` pair via
`NetworkState`. A component never sees the merged id, only its own raw one.

Every node carries:

- **`(P, h)`** — pressure `[Pa]` and specific enthalpy `[J/kg]`, the two
  state variables every fluid property call is built from.
- **`mdot`** — mass flow rate `[kg/s]`, tracked per node independently of
  `(P, h)`.
- **a fluid** — which `BaseFluid` model is actually flowing through that
  node (see *Composition propagation* below); the overwhelming majority of
  nodes just inherit the network's one default fluid.

### Configuring parameters

Constructor kwargs are fine for the parameters you set at creation time, but
a component's tunable parameters can also be set afterwards — including
after it's already been added to a `Network` — via `component.set(**kwargs)`
(`BaseComponent.set()`, inherited by every component):

```python
turb = Turbine(name="turb", map_path="T100 Turb.tur")
network.add_component(turb)

turb.set(N=65000.0)
```

This is useful once a component has several parameters: setting them one
line at a time (rather than all packed into one constructor call) makes it
easy to scan a script later for exactly what was configured. `set()`
validates keys against attributes that already exist on the instance, so a
typo raises `AttributeError` instead of silently creating a new attribute,
and it auto-invalidates the owning network's caches if the component has
already been added. Plain attribute assignment (`turb.N = 65000.0`) still
works identically — `set()` is a validated, typo-safe spelling of the same
thing, not a replacement mechanism.

### Residual equations

`residuals(state) -> list[float]` is the one method every component must
implement. Each entry is written as `computed_value - target_value`, zero
at the solution. A typical two-port flow component contributes exactly
three: a **momentum** residual (the pressure relationship across it), an
**energy** residual (the enthalpy relationship), and a **mass** residual
(`mdot_out == mdot_in`, unless the component itself adds or removes mass,
like a combustor's fuel).

`Network.solve()` assembles one flat vector of unknowns — `(P, h)` for every
node not fixed by a boundary component, `mdot` for every node not fixed, and
one unknown per free parameter (below) — and one flat vector of residuals,
the concatenation of every component's `residuals()`. Newton-Raphson needs a
**square** system; a mismatch raises `NetworkTopologyError` naming the
imbalance, and `Network.check_wiring()` can name the *specific* free
parameter nothing closes (or the two components fighting over one) before
you even call `solve()`.

### Free parameters: giving a target instead of a direct input

Many components accept `None` for one of their own constructor arguments —
a `Compressor`'s shaft speed `N`, a `Combustor`'s `mdot_fuel` — to mean "solve
for this instead of fixing it." That turns the quantity into an extra Newton
unknown (`free_parameters() -> dict[str, guess]`) that needs **exactly one**
matching residual from somewhere else in the network. That residual almost
always comes from one of the components in [Control &
instrumentation](control/index.md):

- **`Setpoint`** ties a free parameter to the *same* component's own
  `report_metrics()` (e.g. a compressor's free `N` to a target power).
- **`Controller`** ties a free parameter to an independent `Sensor`
  reading elsewhere in the network (e.g. a combustor's free `mdot_fuel` to
  a downstream temperature sensor) — an ideal, infinite-gain controller
  with no dynamics, appropriate for a steady operating point.
- **`PIDController`** is the finite-response, time-domain counterpart, for
  use inside `Network.solve_transient()`.
- **`Shaft`** closes the speed-tie residuals between two or more
  turbomachines sharing a physical shaft.

### Differential parameters: state that evolves in time

A handful of components own genuine time-domain state instead of leaning on
an external closing residual: a dynamic `Shaft`'s rotor speed, a `Tank` or
`Drum`'s pressure/enthalpy, a `ThermalMass`'s temperature.
`differential_parameters() -> dict[str, initial]` declares the unknown (same
mechanism as `free_parameters()`), and `state_derivative(state) -> dict[str,
rate]` gives its rate of change. The solver closes it automatically:

- **`Network.solve()`** (steady state): `state_derivative() == 0` — the
  equilibrium where nothing is actually changing.
- **`Network.solve_transient()`**: backward-Euler, `(value -
  value_at_previous_step) / dt == state_derivative()`.

### Reporting

`report_metrics(state) -> dict[str, float] | None` exposes named,
unit-labeled results (`"power [W]"`, `"PR [-]"`, ...) for
`SolveResult.print_report()`, grouped by `report_category()` (e.g.
`"turbomachinery"`, `"heat_exchanger"`, `"controller"`). Components that
return `None` from both (e.g. `Source`, `Sink`, `Pipe`) are left out of the
categorized tables entirely. This is also how `Setpoint`/`Controller` read
the quantity they're driving — through the *same* `report_metrics()` dict
every printed table uses, so there's no separate reporting-vs-control code
path to keep in sync.

### Units

Every component and the solver itself work in strict SI (`Pa`, `K`, `J/kg`,
`kg/s`, `W`). `Source`/`Sink` are the exception: their `P`/`T` constructor
arguments are converted from whatever unit the module-level `settings`
singleton (`thermowave.core.settings`) is configured for (`Pa`/`kPa`/`MPa`/
`bar`/`atm` for pressure, `K`/`C` for temperature) — everything else in the
network downstream of them is already SI.

### Composition propagation

A `Network` has one shared `BaseFluid` by default, but `NetworkState.fluid_at(node)`
lets a component's outlet carry a *different* fluid than its inlet — the
mechanism [`Combustor`](combustion/combustor.md) uses to feed real
equilibrium-product composition downstream, and [`Junction`](flow-elements/junction.md)
uses to blend two differently-composed inlets. Every component that reads fluid properties
reads through `fluid_at()`, so a composition change anywhere upstream is
visible to everything downstream automatically. See each component's own
page for where this matters.

## Diagram legend

Every component page below shows a labeled P&ID-style diagram using this
convention:

- **Solid black arrows** — flow ports (`in`/`out`, or named streams like
  `hot_in`/`cold_out`), the shared `(P, h, mdot)` nodes `connect()` merges.
- **Dashed red arrows** — a `heat_path` (`Convection`/`Conduction`/
  `Radiation`), an optional coupling to the [thermal network](thermal-network/index.md).
- **Blue stub with a dot** — a mechanical shaft connection (`N`), tied to
  other machines via [`Shaft`](mechanical/shaft.md).
- **Dashed rounded badge** — a free parameter (`N`, `mdot_fuel`, ...) left
  `None` at construction, waiting for a `Setpoint`/`Controller`/`Shaft`
  elsewhere in the network to pin it down.
