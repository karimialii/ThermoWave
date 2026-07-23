# ThermoWave

A headless, 1D implicit thermodynamic network solver for steady-state
thermal-fluid systems. You build a network out of components (sources, pipes,
valves, compressors, turbines, heat exchangers, combustors, shafts, sensors,
controllers, ...), wire their ports together, and call `solve()`. Internally
every component contributes unknowns and residual equations to one square
system, solved by damped Newton-Raphson with a finite-difference Jacobian.

```python
from thermowave.fluids.ideal_gas import IdealGasFluid
from thermowave.core.network import Network
from thermowave.components.source import Source
from thermowave.components.pipe import Pipe
from thermowave.components.sink import Sink

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
net = Network(fluid=air)

src = Source(name="src", P=200000.0, T=300.0, mdot=1.0)
pipe = Pipe(name="pipe", L=10.0, D=0.1, f=0.02)
sink = Sink(name="sink")

net.add_component(src)
net.add_component(pipe)
net.add_component(sink)
net.connect(src, "out", pipe, "in")
net.connect(pipe, "out", sink, "in")

result = net.solve()
result.print_report()
```

The full network model — every component, the solver, transient simulation,
plotting, and the two-phase/Rankine machinery — is covered section by
section below. See `tests/` for runnable, working usage of every component
and network pattern described here.

## Install

```bash
pip install -e .
pip install -e ".[coolprop]"   # real-fluid properties via CoolProp
pip install -e ".[cantera]"    # equilibrium-chemistry Combustor + CanteraFluid
pip install -e ".[full]"       # everything, including dev/test tooling
```

## Architecture

### Fluids

`BaseFluid` (`src/thermowave/fluids/base_fluid.py`) is the property interface
every component calls into: `enthalpy_pt(P, T)`, `temperature_ph(P, h)`,
`cp(P, T)`, etc. Four implementations exist today:

- `IdealGasFluid` (`fluids/ideal_gas.py`) — constant `R`/`cp` for a single
  named substance, closed-form properties, no external dependency.
- `IdealGasMixtureFluid` (`fluids/ideal_gas_mixture.py`) — the same
  constant-cp simplification, but blended from named species by mass
  fraction (`{"N2": 0.767, "O2": 0.233}`) instead of one fixed `(R, cp)`
  pair — humid air, flue gas, a fuel-air mix, whatever composition you name,
  with no extra dependency. `IdealGasMixtureFluid.SPECIES` has a small
  built-in table (molar mass + constant molar cp); pass `extra_species=` to
  add/override entries.
- `CoolPropFluid` (`fluids/real_fluid.py`) — wraps
  [CoolProp](http://www.coolprop.org/) for real single-species fluid
  properties (temperature-dependent, not constant-cp). Requires the
  `coolprop` extra.
- `CanteraFluid` (`fluids/cantera_fluid.py`) — real (temperature-dependent)
  thermo for an arbitrary named *mixture*, backed by
  [Cantera](https://cantera.org/) — `IdealGasMixtureFluid`'s use case, but
  without the constant-cp approximation. Requires the `cantera` extra.

A `Network` has exactly **one** fluid instance shared by every component in
it — there's no multi-fluid mixing or composition change *within* a network
today (a `Combustor`'s exhaust, for instance, still reports its outlet
temperature back onto the network's one `BaseFluid`, not a genuinely
different combustion-product mixture — see `Combustor.product_composition()`
below for how to see that chemistry anyway, just not fed back into the
network state).

### Network graph

`Network` (`src/thermowave/core/network.py`) holds the component list and a
union-find over port ids. `add_component()` registers a component's ports
(and any `internal_nodes()`) as graph nodes; `connect(comp_a, "out",
comp_b, "in")` merges two components' local port ids into one shared
network node — after that, both components read/write the same `(P, h)` in
`NetworkState`. Only `kind="flow"` connections exist today (shared
pressure/enthalpy node); `kind="mechanical"` / `"signal"` / `"heat"` are
reserved names but not implemented — instead, cross-component coupling that
doesn't fit the flow-node model (shaft speed matching, controller setpoints)
is expressed as an ordinary component with its own `residuals()` that reads
other components' `report_metrics()` or `free_parameters()`. `Shaft`,
`Controller`, and `Setpoint` are the existing examples of this pattern — see
"Cross-component coupling" below before reaching for a new `Connection` kind.

### The solver

`Solver.solve()` (`src/thermowave/core/solver.py`) assembles one flat vector
of unknowns:

1. `(P, h)` for every network node not fixed by a boundary component (2
   unknowns/node).
2. `mdot` for every node not fixed by a boundary component's
   `fixed_node_mdot()`.
3. One unknown per `(component, key)` pair returned by any component's
   `free_parameters()`.

and one flat vector of residuals: the concatenation of every component's
`residuals(state)` list, in component-add order. Newton-Raphson requires a
**square** system — equal unknown and residual counts — or `solve()` raises
`NetworkTopologyError` before attempting to solve. This is almost always the
error you'll hit while wiring up a new network: it means some free parameter
has no matching residual (or vice versa) somewhere in the graph.

Before solving, the solver forward-propagates a rough `(P, h)` guess through
the network via each component's `guess_outlet()`, and lets each component
refine its own free-parameter guesses via `guess_free_parameters()` using
that warm-started inlet state. This exists because a flat default guess
(every free node starting at the network's first fixed boundary state) can
put a map-based component so far outside its valid operating range at
iteration 0 that the very first Jacobian is singular — see the long comment
in `Solver.solve()` for the full reasoning.

### Components

Every component implements `BaseComponent`
(`src/thermowave/components/base_component.py`). See "Writing a new
component" below for the full contract.

### Reporting

`SolveResult.print_report()` (`src/thermowave/core/reporting.py`) prints one
table per component category (declared via `report_category()`) followed by
a per-node `(P, T, h, mdot)` table. Categories in display order today:
`turbomachinery`, `heat_exchanger`, `combustor`, `controller`, `shaft`,
`sensor`, `generator`. A component with `report_category() -> None` (the
default — `Source`/`Sink`/`Pipe`/`Valve`/`Junction`) is left out of the
categorized tables entirely.

### Units and shared constants

`src/thermowave/core/settings.py` holds the module-level `settings`
singleton (`pressure_unit`, `temperature_unit`) that `Source`/`Sink`
consult when converting constructor inputs to SI — every other component
and the solver itself work in strict SI regardless. Pressure accepts
`Pa`/`kPa`/`MPa`/`bar`/`atm`; temperature only `K`/`C` — Fahrenheit isn't
supported, deliberately (see `temperature_to_si()`'s docstring), since
nothing downstream needs it and it's an easy source of mixed-unit bugs to
add without a concrete use case.

`src/thermowave/core/constants.py` holds physical/unit constants shared by
more than one module (e.g. standard atmosphere, `Pa` per `bar`) — added so
a value used in two or more places has exactly one definition instead of
being copy-pasted (and able to silently drift) across each call site.

## Writing a new component

A component is any class implementing `BaseComponent`. The two required
methods:

- **`ports() -> dict[str, str]`** — named ports (`"in"`, `"out"`, ...)
  mapped to this component's own local port ids, conventionally
  `f"{self.name}.{port_name}"`. `Network.connect()` merges two components'
  port ids into one shared node; a component never sees the merged id
  itself, only its own raw one, when it calls `state.node(self._inlet_node)`.
- **`residuals(state) -> list[float]`** — the equations this component
  contributes, each written as `computed_value - target_value` (zero at the
  solution). Read node state via `state.node(name) -> (P, h)`,
  `state.mdot(name) -> float`, and `state.fluid` for property calls. A
  simple two-port component typically contributes 3 residuals: momentum
  (pressure relationship), energy (enthalpy relationship), and mass
  (`mdot_out == mdot_in`).

Everything else on `BaseComponent` is optional and defaults to "this
component doesn't need it":

| Hook | Default | Override when... |
|---|---|---|
| `fixed_node_values(fluid)` | `{}` | this component is a boundary condition fixing `(P, h)` on one of its ports (e.g. `Source`). |
| `fixed_node_mdot()` | `{}` | this component fixes `mdot` on one of its ports (e.g. `Source`). |
| `internal_nodes()` | `[]` | this component needs solver-tracked nodes that aren't exposed as ports (e.g. a multi-element `Pipe`'s interior discretization points). |
| `free_parameters() -> dict[str, guess]` | `{}` | this component has a scalar unknown of its own, not tied to any node (e.g. a `Compressor`'s shaft speed `N` when left `None` to be solved for instead of given directly). Add exactly one matching residual in `residuals()`, read back via `state.param(f"{self.name}.<key>")`. Keep the guess here cheap/generic — a hardcoded reference value is fine. |
| `guess_free_parameters(fluid, P_in, h_in, mdot)` | falls back to `free_parameters()` | a better free-parameter guess is possible using this component's actual (warm-started) inlet state — e.g. deriving a shaft-speed guess from the real inlet temperature rather than an assumed reference temperature. |
| `guess_outlet(P_in, h_in, mdot)` | pass through unchanged | this component has a large, predictable `(P, h)` swing across it (a compressor's pressure ratio, a heat exchanger's duty) — without this, downstream free nodes can start Newton so far from the solution that the first Jacobian is singular. |
| `report_metrics(state) -> dict \| None` | `None` | this component should show up in a report table — return a `{"name [unit]": value}` dict. |
| `report_category() -> str \| None` | `None` | pick an existing category (see `reporting._CATEGORY_TABLES`) or a new one; components sharing a category share one column layout, so keep `report_metrics()` keys consistent within a category. |
| `differential_parameters() -> dict[str, initial]` | `{}` | this component owns a scalar quantity that evolves in time (e.g. `Shaft(dynamic=True)`'s rotor speed). Like `free_parameters()`, but closed automatically by the solver (see "Transient simulation" below) instead of needing a matching residual — don't add one yourself. |
| `state_derivative(state) -> dict[str, rate]` | `{}` | required alongside `differential_parameters()`: `d(value)/dt` for each declared key, evaluated at the current state. |

### Cross-component coupling (no new `Connection` kind needed)

Some components don't move fluid at all — they exist purely to tie other
components' free parameters or metrics together. `Shaft`
(`src/thermowave/components/shaft.py`) ties two or more components' shaft
speeds together via residuals of the form `component[i].N - gear_ratio *
component[0].N`; `Controller`/`Setpoint` tie a component's free parameter to
a `Sensor` reading or the component's own metric. These have `ports() ->
{}` (or ports with no flow significance) and no `fixed_node_*`; they exist
in `network.components` purely to add unknowns/residuals or read other
components' `report_metrics()`. This is the established pattern for anything
that isn't a flow connection — reach for it before adding a new
`Connection` kind.

### A minimal example

```python
class Valve(BaseComponent):
    def __init__(self, name: str, K: float):
        self.name = name
        self.K = K
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def residuals(self, state: "NetworkState") -> list[float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        mdot = state.mdot(self._inlet_node)

        momentum_residual = (P_in - P_out) - self.K * mdot**2
        energy_residual = h_out - h_in  # isenthalpic
        mass_residual = state.mdot(self._outlet_node) - mdot
        return [momentum_residual, energy_residual, mass_residual]
```

## Testing

```bash
pytest -q
```

212+ tests, with CoolProp/Cantera-gated tests auto-skipped
(`pytest.importorskip`) when those optional dependencies aren't installed.

## Transient simulation

Steady state and transient are the same Newton system — the only thing that
changes is how a *differential* unknown's closing equation is supplied.
Any component can own a scalar piece of time-domain state by overriding two
`BaseComponent` hooks:

- **`differential_parameters() -> dict[str, float]`** — `param_name ->
  initial value`, one Newton unknown per entry (same mechanism as
  `free_parameters()`).
- **`state_derivative(state) -> dict[str, float]`** — `d(value)/dt` for each
  of those keys, evaluated at the current state.

Unlike `free_parameters()`, a differential parameter doesn't need an
external residual to close it (no `Setpoint`/`Controller` required) — the
solver supplies it automatically:

- **`Network.solve()`** (steady state): `state_derivative() == 0` — i.e.
  steady state is the equilibrium where nothing is actually changing.
- **`Network.solve_transient(duration, dt, initial=None, ...)`**
  (`src/thermowave/core/transient.py`): backward-Euler, `(value -
  value_at_previous_step) / dt == state_derivative()`, one Newton solve per
  timestep. `initial` is a `SolveResult` to start from — a prior
  `Network.solve()` at a known operating point, the last step of an earlier
  `solve_transient()` run (to continue it), or omitted entirely, in which
  case `solve_transient()` runs an ordinary steady-state `Network.solve()`
  first to get a genuine equilibrium initial condition.

**Adaptive time-stepping**: `solve_transient(..., adaptive=True, rtol=1e-3,
atol=1e-6, dt_min=None, dt_max=None, safety=0.9, growth_limit=5.0,
shrink_limit=0.2, max_step_shrinks=10)` — `dt` becomes the *initial* step
size rather than a fixed one. Each step is tried once at the current size
`h` and once as two `h/2` half-steps (step-doubling); since backward-Euler's
local error is O(dt²), the difference between the two estimates the
half-step (more accurate) result's own error. That result is kept if the
weighted-RMS error across every differential state is within tolerance
(`rtol`/`atol` combine the same way as `scipy.integrate.solve_ivp`: `scale =
atol + rtol*|value|`), and the next step's size is rescaled accordingly
(`safety*(1/err)**0.5`, clamped to `[shrink_limit, growth_limit]*h` and to
`[dt_min, dt_max]`). A rejected step (error over tolerance, or the Newton
solve itself failing to converge at that `h`) retries smaller without
advancing time; `max_step_shrinks` caps retries before raising
`ConvergenceError` rather than looping forever on a genuinely unsolvable
step. Costs roughly 3x the nonlinear solves of a fixed-step run of the same
duration — worth it for a network whose time constants vary enough across
the run (e.g. a fast initial transient settling into a slow one) that no
single fixed `dt` is efficient throughout. Requires at least one component
declaring `differential_parameters()` (step-doubling's error estimate needs
something to compare) — a network driven only by `step()`-able components
(e.g. a lone `PIDController`) must use `adaptive=False`. See
`tests/test_transient.py` for worked examples, including the `dt_max` clamp
and the `step()`-called-once-per-accepted-step behavior with
multi-component networks.

`solve_transient()` takes no `shaft`/`setpoint`/similar arguments — it
discovers what to integrate purely from what the network's own components
declare, so it's the same code path for any current or future component
with its own time constant. It raises `ValueError` if nothing in the
network has time-varying state at all (differential parameters or a
`step()`-able component — see `PIDController` below): there'd be nothing
for it to do that a plain `Network.solve()` doesn't already cover.

Only declared differential state (and anything a `step()`-able controller
drives) carries a real time history — there's no fluid inertia/accumulation
term, so every step is otherwise a genuine steady state. This holds well
when the physical time constant you're modeling (rotor inertia, a
controller's actuator) is much slower than the flow/thermal ones, which is
typical for gas-turbine spool-up/load-step and control-loop transients.

**Rotor spin-up**: `Shaft(..., dynamic=True, inertia=..., N0=...)` (`src/
thermowave/components/shaft.py`) owns its own speed as differential state,
integrated from net torque: `d(N)/dt = (net_power / omega) / inertia`.
Every connected component's own free speed parameter (e.g.
`Compressor(..., N=None)`) is tied to the shaft's speed by `Shaft`'s own
residuals — no `Setpoint` needed anywhere. `dynamic=False` (the default)
keeps `Shaft`'s original steady-state-only behavior: components' speeds are
tied to `components[0]`'s, which is a plain algebraic free parameter still
needing an external `Setpoint`/`Controller` to pin down — appropriate for a
governed engine, where speed is actively regulated rather than settling
under its own inertia. See `tests/test_transient.py` for a dynamic
spin-up example and `tests/test_shaft.py` for the governed case.

**Volume dynamics**: `Tank(name, V, P0, T0, fluid, heat_loss=None)`
(`src/thermowave/components/tank.py`) owns its own contents' `(P, h)` as
differential state — the genuinely transient counterpart to `Junction`'s
zero-volume, quasi-steady mixing. Everything else in the network is an
algebraic function of its inlet at the current instant; a `Tank` actually
accumulates, so `mdot_in` and `mdot_out` can (and generally do) differ while
its own pressure/temperature lag behind — real filling, blow-down, and surge
behavior instead of only ever a sequence of independent steady states. Mass
and energy conservation on the contents are solved as a 2x2 linear system
for `dP/dt`, `dh/dt`, using finite-differenced `fluid.density_ph()` partials
(the same finite-difference philosophy the solver's own Jacobian uses)
rather than assuming a closed-form ideal-gas equation of state — works with
any `BaseFluid`. Neither `mdot` is constrained by the tank itself (a real
control volume's in/outflow are set by whatever's upstream/downstream, not
the volume); pair its outlet with a `Valve` or similar flow-resistance
component to close the system, the same way a `Source(mdot=None)` needs
something downstream to pin its own flow.

**Closed-loop control**: `PIDController`
(`src/thermowave/components/pid_controller.py`) is the time-domain
counterpart to `Controller`/`Setpoint`, which pin a target exactly on every
steady-state solve (an ideal, infinite-gain controller with no dynamics).
`PIDController` pins its actuated free parameter to `self.output` on every
solve (same residual shape as `Controller`), but `self.output` only
*changes* once per transient step, via `PIDController.step()` computing
`output0 + Kp*error + Ki*integral + Kd*derivative` from the current
`Sensor` reading — so it takes several timesteps to settle, with real
overshoot/offset/settling-time behavior. `step()`-able components (found
by simple duck-typing, `hasattr(component, "step")`) are discovered and
stepped automatically by `solve_transient()`, same as differential
parameters — no separate argument needed. See `tests/test_pid_controller.py`
for single- and multi-loop usage.

## Plotting

`thermowave.core.plotting.ThermoPlot` (requires the optional `plot` extra —
`pip install thermowave[plot]`, matplotlib) is a themed, chainable plotting
class: a colorblind-safe palette, light background, subtle gridlines, and
trimmed spines applied per-figure (never mutates matplotlib's global
rcParams), with methods for every chart type this kind of analysis needs —
`.line()`, `.scatter()`, `.series()` (component-metric sweeps), `.twin_axis()`
(two independent y-axes sharing one x-axis), `.transient()` (time series),
`.map()` (a `CharacteristicMap`'s iso-speed curves with an operating-point
overlay), and `.bar()`. Mutating calls return `self`, so they chain:

```python
from thermowave.core.plotting import ThermoPlot

ThermoPlot(title="Shaft spool-up", ylabel="N [rev/min]").transient(
    history, [(shaft, "N [rev/min]")]
).finish()

power = ThermoPlot(xlabel="mdot [kg/s]", ylabel="power [W]")
power.line(mdots, powers, label="power [W]")
power.twin_axis(ylabel="eta_s [-]").line(mdots, etas, label="eta_s [-]")
power.finish()  # one combined legend across both axes

ThermoPlot(title="Compressor map").map(
    comp.map, kind="pressure_ratio", operating_points=[(B, PR), ...]
).finish()
```

`.ax`/`.fig` are public attributes for any matplotlib call the class
doesn't wrap directly. See `tests/test_plotting.py` for usage of each
chart type.

`TransientResult.plot(*series, ...)` and `plotting.plot_series(x, results,
series, ...)` still work exactly as before — they're now thin wrappers
around `ThermoPlot.transient()`/`.series()` respectively, so every existing
call site gets the new theme automatically:

```python
history = network.solve_transient(duration=1.0, dt=0.05)
history.plot((shaft, "N [rev/min]"), title="Shaft spool-up")
history.plot((comp, "power [W]"), (turb, "power [W]"), ylabel="power [W]")
```

`plot_series()` is generic over any list of `SolveResult` with a matching
x-axis, so it also covers manual parameter sweeps — pass your own `x`
values and the `SolveResult` from each `network.solve()` call. Pass `ax=`
to draw onto an existing `Axes` (for subplots), or `show=False,
save_path=...` to save headlessly instead of popping up a window.

## Heat transfer

`src/thermowave/components/heat_transfer.py` adds real thermal-mass and
heat-path physics on top of what was previously fully adiabatic
turbomachinery/combustion (`Pipe`'s own fixed `heat_loss` was the only
exception): `ThermalMass` is a solid's own temperature as a differential
state (a casing, a shaft, ...); `Convection`, `Conduction`, and `Radiation`
are heat paths between two temperature sources. Every `Q(state)` a path
returns is positive when its `a` endpoint is hotter than `b` (heat flowing
`a -> b`), the same sign convention `Pipe`'s own `heat_loss` already uses.

- `ThermalMass(name, thermal_capacitance, T0)` — `thermal_capacitance` is
  the lumped `m*cp` [J/K], given directly (like `MultiPassHeatExchanger`'s
  `UA`, not derived from separate mass/material-cp values). No flow ports —
  it only participates in the differential-state/Newton bookkeeping, the
  same mechanism `Tank`/`Drum`/`Shaft(dynamic=True)` already use.
  `Network.solve()` closes its temperature to whatever value makes its net
  heat zero (steady state); `Network.solve_transient()` integrates it
  forward from `T0` instead.
- `Convection(name, a, b, h, A)` — `Q = h*A*(T_a - T_b)`, covering both
  free and forced convection (they differ only in how `h` is physically
  obtained, not in this formula).
- `Conduction(name, a, b, k, A, L)` — `Q = (k*A/L)*(T_a - T_b)`, steady 1D
  conduction (e.g. through a shaft between two casings).
- `Radiation(name, a, b, emissivity, A, view_factor=1.0)` — `Q =
  emissivity*view_factor*sigma*A*(T_a^4 - T_b^4)`. General-purpose
  surface-to-surface (or surface-to-ambient) physics, not combustor-
  specific — the intended primitive a future 1D combustion-chamber liner
  model will discretize into many of, rather than needing its own formula.

`a`/`b` on any path is one of: a `ThermalMass`, a fixed `float` (e.g.
ambient temperature), or a `(component, port_name)` tuple reading that
component's live fluid temperature at that port's node.

`Turbine`, `Compressor`, `SimpleTurbine`, `SimpleCompressor`,
`SimpleCombustor`, and `Combustor` each have an optional `heat_path`
attribute (`None` by default — fully adiabatic, unchanged from before this
existed) that, when set, actually perturbs that component's own energy
residual — real two-way coupling, not just reporting. Since a path needs
`(component, "out")` as one of its own endpoints, it can only be built
*after* that component already exists, and the component needs to know
about the path too (for its residual) — so wiring one in is a two-step,
plain-attribute-assignment process, the same pattern `Shaft`/`Generator`
already use for referencing components built earlier, extended one step
further because both endpoints need each other here:

```python
turb = Turbine(name="turb", map_path="T100 Turb.tur", N=65000.0)
casing = ThermalMass(name="turb_casing", thermal_capacitance=200.0, T0=300.0)
conv = Convection(name="turb_conv", a=(turb, "out"), b=casing, h=50.0, A=0.3)
to_ambient = Convection(name="turb_ambient", a=casing, b=288.15, h=10.0, A=1.0)

turb.heat_path = conv                                    # turbine's fluid loses Q
casing.heat_sources = [(conv, 1.0), (to_ambient, -1.0)]   # casing gains conv's Q, loses to_ambient's
```

See `tests/test_heat_transfer_integration.py` for the full picture: a
turboshaft where each machine convects into its own casing, the two
casings conduct through a shaft `ThermalMass` between them, both casings
convect to ambient, and the casings' temperatures actually change the
cycle's predicted `T_out`.

## Seeing combustion products

`Combustor.product_composition(state)` (`src/thermowave/components/
combustor.py`) re-runs the same Cantera equilibrium calculation `residuals()`
already uses to find `T_out`, but returns the full result — a `{species:
mole_fraction}` dict — instead of discarding everything but temperature. The
major species (CO2, H2O, O2, N2, plus CO/NO if the mechanism produces them
above a trace threshold) are also surfaced directly in `report_metrics()` as
`"X_<species> [-]"`, so they show up in the printed Combustors table without
calling anything extra. See `tests/test_combustor.py` for usage.

This is genuine product chemistry (dissociation, excess-air/equivalence-ratio
effects, trace NO from GRI-Mech, ...), and — when `Combustor`'s own inlet
fluid is a `CanteraFluid` — it now *does* feed back into the network's node
state: every node downstream of the combustor carries the actual reacted
product composition, not a copy of the pre-combustion working fluid. See
"Composition-aware fluid propagation" below for how that works and its one
real restriction (the inlet fluid has to already be Cantera-based).

## Composition-aware fluid propagation

`NetworkState` now carries an optional per-node fluid map alongside its
per-node (P, h): `state.fluid_at(node)` returns whatever `BaseFluid`
actually reached that node, falling back to the network's own default
`fluid` for every node nothing changed (the overwhelming majority — this
costs nothing for a network with no composition-changing component). Every
component that reads fluid properties (`Pipe`, `Valve`, `Compressor`,
`Turbine`, both heat exchangers, `Tank`, `Nozzle`, `CheckValve`, `Sensor`,
...) reads through `fluid_at()` instead of the network's single `fluid`
directly, so a composition change anywhere upstream is automatically
visible to everything downstream's actual physics, not just its reported
metrics.

`Combustor` is the one component that changes composition today: its
`outlet_fluid()` hook returns the real Cantera equilibrium-product mixture
(the same one `residuals()` already computes for `T_out`, cached per Newton
residual evaluation so both share one `equilibrate()` call) whenever its own
inlet fluid is a `CanteraFluid` — one call to
`Network._resolve_node_fluid()` per residual evaluation forward-propagates
that through the rest of the network via a fixed-point loop over each
component's `warm_start_pairs()`/`outlet_fluid()`, the same order-
independent technique `Solver.solve()`'s own (P, h) warm-start guess already
uses. This stays correct through every Newton iteration even though the
product composition itself depends on a free unknown (`mdot_fuel`).

The one real restriction: this is only physically consistent when the
combustor's inlet fluid already shares Cantera's absolute (formation-
enthalpy-referenced) datum. Mixing that with e.g. `IdealGasFluid`'s
`h = cp*T` (referenced to h=0 at T=0K, with no notion of chemical energy at
all) would silently corrupt any downstream energy balance — so for any
inlet fluid that isn't a `CanteraFluid`, `outlet_fluid()` returns `None`
(pass-through) and `Combustor` falls back to its original behavior:
chemistry-informed `T_out`, tracked downstream through the same fluid model
as upstream.

`Junction` does real composition mixing too, via a second hook
(`merge_fluids()`, alongside `outlet_fluid()`) built specifically for
multi-inlet merge points that need every inlet's resolved fluid at once
rather than one at a time: if every inlet is already the same fluid object
(the common case), it just passes through; if inlets genuinely differ and
each one exposes a `mass_fractions()` method and a matching `mechanism`
(both `CanteraFluid` and `Combustor`'s own equilibrium-product fluid do —
checked structurally, not via `isinstance(..., CanteraFluid)`, precisely so
two *different* combustor exhaust streams merging — the realistic case this
exists for — actually qualifies), it computes a genuine mass-weighted blend
of their compositions. Anything else (fluids that don't expose that
contract, or share it on different mechanisms) falls back to the first
inlet's fluid — an explicit, documented simplification, not a silent
guess. A recycle/EGR loop (composition depending on itself through a flow
cycle) still isn't handled — the fixed-point loop only resolves a node once
its inlet is already known, so a true cycle just never settles that node
(it silently stays at the pass-through default instead); out of scope for
now, unreachable by any example/test network here since nothing recycles.

## Two-phase components (evaporator / condenser / drum / Rankine)

Phase-change equipment is modeled via **saturation enthalpies and quality**,
never `cp` — the effectiveness-NTU `C = mdot*cp` framework the gas heat
exchangers use breaks down in the two-phase dome, where `cp -> infinity`.
This needs a real-fluid saturation model, so it is gated on a
two-phase-capable fluid (`CoolPropFluid`, `pip install thermowave[coolprop]`).

**Fluid foundation.** `CoolPropFluid` gained saturation/quality methods
(`saturation_temperature`, `saturated_liquid_enthalpy`/`saturated_vapor_enthalpy`,
`enthalpy_pq`, `quality_ph`, `saturation_pressure`) and entropy methods
(`entropy_ph`, `enthalpy_ps`). These are **additive to `CoolPropFluid` only**,
never to the `BaseFluid` interface (the ideal-gas/Cantera models can't
implement them). Components detect the capability structurally, via
`thermowave.fluids.two_phase.supports_two_phase()` / `require_two_phase()`
(and `require_entropy()`), the same duck-typing rationale `Junction.merge_fluids`
uses — a future REFPROP-backed fluid could satisfy the contract without
subclassing.

**Components** (all take the phase-change side's duty from the specified
outlet state, e.g. saturated vapor / a quality / N degrees of superheat, then
`Q = mdot*Δh`):

- **`SimpleEvaporator` / `SimpleCondenser`** — single-stream (boil/condense to
  a spec outlet, or a given duty), the SimpleCombustor-style model that
  doesn't model the heat source/sink explicitly. 3 residuals.
- **`Evaporator` / `Condenser`** — two-stream, coupling the phase-change side
  to an explicit heat-source / coolant stream (SimpleHeatExchanger-style, 4
  ports / 6 residuals). Both streams share the network's single fluid model
  (same limitation as `SimpleHeatExchanger`), so these model a same-fluid
  coupling. The stream-to-stream pinch is a **reported diagnostic, not a
  solved constraint**: a negative reported `pinch [K]` flags an infeasible
  outlet spec (the source would have to end up colder than the boiling fluid)
  rather than being silently corrected — the user owns feasibility, the way
  `SimpleHeatExchanger` trusts a user-supplied effectiveness.
- **`Drum`** — a steam drum: the two-phase, multi-port analogue of `Tank`
  (differential `(P, h)` on a saturated inventory via the identical
  finite-difference 2×2 mass/energy ODE, generalized to feed/riser inlets and
  saturated-vapor/saturated-liquid outlets). A drum's level is a **pure
  integrator** with no steady-state restoring force — real drums need level
  control, and correspondingly a plain steady `Network.solve()` is *singular*
  in the drum level (there is no algebraic level to solve for). The Drum is
  therefore a transient component: its `state_derivative()` is what
  `Network.solve_transient()` integrates. See `tests/test_drum.py` for the
  level/pressure response to a steam-demand step (including the pressure
  "shrink" on subcooled-feed injection).
- **`Pump` / `SteamTurbine`** — the entropy-based isentropic-path components a
  closed steam cycle needs. `Pump` closes the low-pressure side back to boiler
  pressure; `SteamTurbine` is the wet-steam-correct counterpart to
  `SimpleTurbine` (whose ideal-gas gamma-relation is wrong once the expansion
  crosses into the dome), and reports exhaust quality `x_out` (blade-erosion
  concern). Both use `entropy_ph`/`enthalpy_ps`.

A full Pump → boiler (superheat) → SteamTurbine → condenser cycle on water
(reporting cycle efficiency and wet-steam turbine exhaust) can be built the
same way as the other network examples above. It is an **open** chain (a
Source pins the feedwater, a Sink terminates the exhaust) — the unrolled
equivalent of the closed cycle, since the solver requires a fixed boundary
node and true recycle loops are out of scope (see the fluid-propagation
section above).

## Roadmap

- **Adaptive time-stepping — landed.** `Network.solve_transient(...,
  adaptive=True)` replaces a fixed `dt` with step-doubling error control
  (`rtol`/`atol`, `dt_min`/`dt_max`, a PI-style step-size controller) — see
  "Transient simulation" above for the full contract. `adaptive=False` (the
  default) is the original fixed-step behavior, unchanged.
- **Two-phase / Rankine components — landed.** `SimpleEvaporator`/
  `SimpleCondenser` (single-stream), `Evaporator`/`Condenser` (two-stream),
  `Drum` (transient steam drum), and `Pump`/`SteamTurbine` — all gated on a
  two-phase-capable `CoolPropFluid` via the duck-typed
  `two_phase.require_two_phase()`/`require_entropy()` helpers. See the
  "Two-phase components" section above for the full model and its scope
  boundaries (pinch is a reported diagnostic, not a solved constraint; a
  drum's level is an integrator so its steady solve is singular by design —
  it's a transient component; two-phase boundary streams can't be pinned
  through a `Source`'s `(P, T)` on the saturation line). Still open: a
  boundary component that fixes a two-phase state directly (P + quality),
  and a UA/pinch-solved two-stream evaporator (rather than duty from the
  outlet spec).
- **Composition-aware fluid propagation — landed**, including real
  `Junction` mixing (a mass-weighted blend of two or more differently-
  composed inlets, gated on all of them exposing a common
  `mass_fractions()`/`mechanism` contract). See its own section above for
  the full mechanism. Still open: recycle/EGR loops — explicitly out of
  scope for now (flagged in that section), unreachable by any example/test
  network today since nothing recycles.
- **More components — landed.** `Tank` (constant-volume plenum with real
  mass/energy storage, see "Volume dynamics" above), `ElectricMotor`
  (electrically-driven mechanical power, the inverse of `SimpleGenerator`),
  `CheckValve` (one-way flow restriction), `Nozzle` (converging nozzle,
  isentropic expansion to velocity, with a choked-flow mass-flow cap), and
  `MultiPassHeatExchanger` (effectiveness derived from UA/flow arrangement
  via the effectiveness-NTU method, instead of taken directly as an input —
  see its own docstring for what "multi-pass" does and doesn't buy you
  here) have all landed. `Nozzle` also now supports an optional diverging
  section (`D_exit`): once choked, flow expands isentropically through it to
  whatever supersonic Mach the geometric area ratio implies (the standard
  area-Mach relation), reported as "Mach_exit [-]"/"P_exit_ideal [Pa]" —
  mdot is still set entirely by throat continuity, unaffected by the
  diverging section, and D_exit is inert while unchoked (see `Nozzle`'s own
  docstring for the full contract and its scope boundary: this is the
  design/shock-free exit condition only, not a shock-train or off-design
  model). `MultiPassHeatExchanger` also gained a real reversing-header
  multi-pass arrangement: `arrangement="shell_and_tube"` treats `n_passes`
  as the number of *shell* passes N (each internally 2 tube passes, the
  standard TEMA "1-2N" configuration) and computes effectiveness via the
  genuine Bowman/Mueller/Nagle F-correction-factor relation (as given in
  Incropera) instead of chaining discretized stages — unlike the other
  three arrangements, where n_passes is mathematically a no-op
  (counterflow/parallel) or not guaranteed to help (crossflow),
  `shell_and_tube`'s n_passes is a real physical lever: effectiveness
  strictly increases with N, converging toward the true counterflow limit
  as N grows, and reduces exactly to a single 1-2 shell-and-tube pass at
  N=1 — see `MultiPassHeatExchanger`'s own docstring for the closed-form
  relations and why it needs no internal discretization nodes (the
  tube-side reversal is captured analytically, not modeled as a node
  chain).

### Future development

Genuine open gaps in the current architecture, not just unwritten
components — each is a real wall a general user hits today, with a pointer
to where it lives:

- **`Network.connect()` only implements `kind="flow"`**
  (`src/thermowave/core/network.py`, `_SUPPORTED_CONNECTION_KINDS`).
  Mechanical coupling (`Shaft`) and heat coupling (`Convection`/
  `Conduction`/`Radiation`) each exist only as their own bespoke component,
  not as a first-class connection kind in the graph API — there's no
  `network.connect(a, "shaft_out", b, "shaft_in", kind="mechanical")`. A
  uniform connection API across domains (flow/mechanical/heat/electrical)
  would let new coupling types plug into the existing graph-traversal and
  reporting machinery instead of each needing its own component class.
- **`CharacteristicMap` only reads single-angle turbomachinery maps**
  (`src/thermowave/maps/characteristic_map.py`) — a variable-geometry
  compressor/turbine map (multiple vane-angle blocks in one `.cop`/`.tur`
  file) isn't parseable at all. Real VGV/VGT hardware needs this; it's a
  file-format gap, not a modeling one — the underlying corrected-flow/PR/eta
  interpolation already generalizes.
- **`MultiPassHeatExchanger` is closed over 4 built-in arrangements**
  (`src/thermowave/components/multi_pass_heat_exchanger.py`,
  `_ARRANGEMENTS`) — `counterflow`/`parallel`/`crossflow`/`shell_and_tube`
  only, enforced by a `ValueError` on anything else. There's no hook for a
  user's own effectiveness-NTU correlation (plate-fin, finned-tube, a
  vendor's empirical curve); a `Callable[[NTU, Cr], float]` override
  alongside the named presets would open this up without disturbing the
  existing four.
- **One fluid per `Network`** (`Network.fluid` in
  `src/thermowave/core/network.py`) — a network has a single default
  working-fluid model, with per-node overrides only reaching as far as
  components that implement `outlet_fluid()`/`merge_fluids()` (combustion
  products, real mixing at a `Junction`). A genuinely independent
  multi-fluid system (e.g. two unrelated streams in a binary/organic-
  Rankine cascade with no shared enthalpy datum) needs separate `Network`
  instances today rather than one connected graph.
- **Recycle/EGR loops are explicitly out of scope** (see "Composition-aware
  fluid propagation" above) — the solver requires an acyclic flow graph
  with a fixed boundary `Source`; a stream whose composition depends on
  itself through a closed loop has no representation yet. This is a bigger
  lift than the others (it likely needs the solver's own unknown-discovery
  pass to handle a genuine fixed-point over composition, not just
  P/h/mdot), so it's listed last.
- **No CI.** Tests and `ruff`/`mypy` (already in the `dev` extra) aren't
  wired into a GitHub Actions workflow yet — right now "does it still pass"
  is only ever checked locally before a push.
- **`pyproject.toml` has no `[project.urls]`** — the PyPI project page has
  no link back to this repository, issue tracker, or docs. Small, but worth
  fixing before the next release rather than accumulating with the rest.
