# ThermoWave Foundation — Design Spec

**Status:** Approved
**Date:** 2026-07-14
**Scope:** First sub-project of the ThermoWave rebuild. Establishes the core architecture — fluid models, network graph, Newton-Raphson solver, and a minimal component set (Source, Sink, Pipe) — sufficient to solve a simple thermofluid network end-to-end. Turbomachinery (compressor, turbine), combustion, heat exchangers, shafts/generators, unit conversion beyond the settings stub, transient simulation, and terminal/plotting reporting are explicitly **out of scope** and will be separate follow-up specs.

## Background

ThermoWave is a from-scratch rebuild of a previously-discarded implementation (same repo, same package name). The prior version's `README.md`/`pyproject.toml` (recovered from git history at commit `b361e8c`, not reused as code) described a 1D implicit thermodynamic network solver for gas turbines and thermal-fluid systems, using CoolProp for real-fluid properties and a staggered node/branch graph solved via `scipy.optimize.root`. This rebuild keeps that overall direction but changes the solver from a scipy black-box to a hand-rolled Newton-Raphson loop for tighter control over damping and physical clamping (CoolProp calls are fragile under wild Newton guesses).

## Goals

- Prove the core architecture (fluid abstraction, graph representation, nonlinear solver, component residual contributions) end-to-end with a minimal but real component set.
- Establish patterns (BaseFluid, BaseComponent, exception types, settings/units) that later specs (turbomachinery, combustion, heat exchangers) will extend without rearchitecting.
- Ship something testable: a `Source → Pipe → Sink` network that solves and matches hand-calculated results.

## Non-Goals

- Turbomachinery, combustion, heat exchangers, shafts, generators, performance maps.
- Transient/time-stepping simulation.
- Multi-fluid networks (mixing, species tracking) — one fluid per `Network` in this spec.
- Rich terminal reporting, plotting, result tables.
- Analytic/sparse Jacobians (finite-difference only for now).

## Architecture

### Package layout

```
src/thermowave/
  __init__.py
  fluids/
    __init__.py
    base_fluid.py       # BaseFluid abstract interface
    ideal_gas.py         # IdealGasFluid
    real_fluid.py         # CoolPropFluid
  core/
    __init__.py
    network.py           # Network
    solver.py            # Solver, Newton-Raphson loop
    settings.py           # settings object, unit conversion helpers
    exceptions.py         # ConvergenceError, FluidRangeError, NetworkTopologyError
  components/
    __init__.py
    base_component.py     # BaseComponent abstract interface
    source.py
    sink.py
    pipe.py
tests/
  test_fluids.py
  test_solver.py
  test_components.py
  test_network_e2e.py
pyproject.toml
README.md
```

### Solver formulation

- Network is a `networkx.DiGraph`. **Nodes** carry scalar thermodynamic state: pressure `P` [Pa] and specific enthalpy `h` [J/kg]. **Branches** (edges) carry mass flow `mdot` [kg/s].
- Each component contributes residual equations evaluated against the current node/branch state (mass balance, momentum balance, energy balance as applicable).
- Boundary conditions (`Source` fixes its node's `P, h`) are **eliminated** from the unknown vector, not added as extra residual rows — the solver never "solves for" a value a `Source` has already fixed. This keeps the system square: `len(unknowns) == len(residuals)`.
- `Solver` assembles the free-variable vector `x` (unfixed node `P,h` pairs + all branch `mdot`), builds `F(x)` from component residual contributions, computes a finite-difference Jacobian `J`, and iterates:

  ```
  x_{k+1} = x_k - damping * solve(J, F(x_k))
  x_{k+1} = clamp(x_{k+1})   # enforce P > P_min, valid fluid range
  ```

  until `‖F(x)‖ < tol` or `max_iter` exceeded (raises `ConvergenceError`).

### Units

- All internal computation is strict SI: Pa, K, J/kg, kg/s, W, m, m².
- `core/settings.py` exposes a module-level `settings` object (e.g. `settings.pressure_unit = 'bar'`) used only for **construction-time convenience and result display** — e.g. `Source(P=1.013, T=15.0)` interpreted in `settings`'s configured units, converted to SI immediately on construction. No component or solver internals ever see non-SI values.

## Fluid Models

`fluids/base_fluid.py` — `BaseFluid` (ABC):
```python
class BaseFluid(ABC):
    name: str

    @abstractmethod
    def density_ph(self, P: float, h: float) -> float: ...

    @abstractmethod
    def temperature_ph(self, P: float, h: float) -> float: ...

    @abstractmethod
    def enthalpy_pt(self, P: float, T: float) -> float: ...

    @abstractmethod
    def cp(self, P: float, T: float) -> float: ...
```

`fluids/ideal_gas.py` — `IdealGasFluid(BaseFluid)`:
- Constructed with `R` (specific gas constant, J/kg·K) and constant `cp` (J/kg·K); `cv = cp - R`.
- `density_ph(P, h)`: derive `T = h / cp` (reference `h=0` at `T=0`), then `ρ = P / (R T)`.
- `temperature_ph(P, h)`: `T = h / cp`.
- `enthalpy_pt(P, T)`: `h = cp * T`.
- No external dependency.

`fluids/real_fluid.py` — `CoolPropFluid(BaseFluid)`:
- Constructed with a CoolProp fluid name string (e.g. `'Air'`, `'Water'`).
- Import of `CoolProp` is deferred to `__init__` — if not installed, raises `ImportError("CoolPropFluid requires the 'coolprop' extra: pip install thermowave[coolprop]")`.
- Wraps `CoolProp.CoolProp.PropsSI('D', 'P', P, 'H', h, name)` etc. for each method.
- Before each `PropsSI` call, clamps `P` and `h` to a configurable valid range (`P_min`, `P_max` defaults e.g. `1e3` to `1e8` Pa) to avoid CoolProp exceptions on wild Newton guesses; if clamping still produces a CoolProp error, re-raises as `FluidRangeError` with the offending values.

## Network & Solver

`core/network.py` — `Network`:
```python
class Network:
    def __init__(self, fluid: BaseFluid): ...
    def add_component(self, component: BaseComponent) -> None: ...
    def solve(self, tol: float = 1e-6, max_iter: int = 100, damping: float = 1.0) -> SolveResult: ...
```
- `add_component` registers the component and wires its declared inlet/outlet node names into the internal `networkx.DiGraph`, creating nodes/edges as needed.
- `solve()` validates topology (raises `NetworkTopologyError` if no `Source` is present, or a node is unreachable from any `Source`), then delegates to `Solver`.

`core/solver.py` — `Solver`:
```python
class Solver:
    def __init__(self, network: Network): ...
    def solve(self, tol: float, max_iter: int, damping: float) -> SolveResult: ...
```
- Builds the free-variable index (which node `P,h` and branch `mdot` values are unknowns vs. fixed by a `Source`).
- `residuals(x)`: scatters `x` back into node/branch state, calls each component's `residuals(state)`, concatenates into `F`.
- `jacobian(x)`: forward-difference, `J[:, j] = (F(x + eps*e_j) - F(x)) / eps`.
- Newton loop with damping and clamping as described above.
- Returns `SolveResult(converged: bool, iterations: int, residual_norm: float, node_states: dict, branch_states: dict)`.

`core/exceptions.py`:
```python
class ConvergenceError(Exception): ...
class FluidRangeError(Exception): ...
class NetworkTopologyError(Exception): ...
```

## Components

`components/base_component.py` — `BaseComponent` (ABC):
```python
class BaseComponent(ABC):
    name: str

    @abstractmethod
    def ports(self) -> dict[str, str]:
        """Named ports -> node names this component connects to."""

    @abstractmethod
    def residuals(self, state: NetworkState) -> list[float]:
        """Residual contributions given current node/branch state. [] if this
        component only fixes boundary values (handled via elimination)."""

    def fixed_values(self) -> dict[str, float]:
        """Node/branch values this component fixes (for elimination). Empty by default."""
        return {}
```

`components/source.py` — `Source(name, P, T)`:
- Single port `outlet`.
- `fixed_values()` returns `{outlet.P: P_si, outlet.h: fluid.enthalpy_pt(P_si, T_si)}` (converted via `settings`).
- `residuals()` returns `[]`.

`components/sink.py` — `Sink(name)`:
- Single port `inlet`.
- No fixed values, no residuals — its presence in the graph closes the branch and satisfies the topology check.

`components/pipe.py` — `Pipe(name, L, D, f, n_elem=1, heat_loss=None)`:
- Ports `inlet`, `outlet`; internally creates `n_elem - 1` intermediate nodes/branches when `n_elem > 1`.
- Per element: momentum residual `P_in - P_out - f*(L/n_elem)/D * (rho*v**2/2) = 0` where `v = mdot / (rho * A)`, `A = pi*D**2/4`, `rho` from `fluid.density_ph` at element-average state.
- Mass residual: `mdot_in - mdot_out = 0` for each element (no storage).
- Energy residual: `h_in - h_out - q_loss/mdot = 0` per element, where `q_loss = 0` if `heat_loss is None`, else `heat_loss` (constant W, split evenly across elements) or `heat_loss(T_element, T_ambient)` if a callable is provided.

## Testing Strategy

- **`test_fluids.py`**: `IdealGasFluid` checked against `PV = mRT` for known `P, T`; `CoolPropFluid` checked against a known reference point (water at 100 °C / 1 atm ⇒ ~2676 kJ/kg saturated vapor enthalpy region — use a safe single-phase point instead, e.g. superheated steam at 200 °C/1 bar), skipped via `pytest.importorskip("CoolProp")` if not installed; clamping behavior tested by passing out-of-range `P`/`h` and asserting `FluidRangeError` or a clamped-but-valid result per the design.
- **`test_solver.py`**: Jacobian correctness on a synthetic 2-equation system with a known analytic Jacobian; convergence on a trivial 1-unknown residual; `ConvergenceError` raised on a deliberately unsolvable system (e.g. `F(x) = 1` for all `x`).
- **`test_components.py`**: `Pipe.residuals` checked against hand-calculated Darcy-Weisbach `ΔP` for one element with fixed `f, L, D, mdot, rho`.
- **`test_network_e2e.py`**:
  - `Source → Pipe(n_elem=1, heat_loss=None) → Sink`: solved `P_out` matches `P_in - ΔP_darcy(hand-calculated)` within `1e-4` relative tolerance; `h_out == h_in` (adiabatic).
  - Same topology with `n_elem=3`: total `ΔP` matches the `n_elem=1` analytic case within tolerance.
  - `Pipe` with constant `heat_loss=Q`: `h_out == h_in - Q/mdot` within tolerance.
  - Topology error case: `Network` with only a `Pipe` and no `Source` raises `NetworkTopologyError` on `solve()`.

## Dependencies (`pyproject.toml`)

```
dependencies = ["numpy>=1.24", "scipy>=1.10", "networkx>=3.0"]

[project.optional-dependencies]
coolprop = ["CoolProp>=6.4"]
dev = ["pytest>=7.4", "pytest-cov>=4.1", "ruff>=0.1", "mypy>=1.5"]
full = ["thermowave[coolprop,dev]"]
```

Note: `scipy` is retained as a dependency for the linear solve step (`scipy.linalg.solve` for `J @ dx = -F`) even though the outer nonlinear loop is hand-rolled, not `scipy.optimize.root`.

## Open Follow-Up Specs (not in this scope)

1. Turbomachinery: `Compressor`, `Turbine` with performance-map integration.
2. Combustion: `Combustor` with stoichiometric/reacting-gas fluid model.
3. Heat exchangers: ε-NTU model, `HeatExchanger` component.
4. Mechanical: `Shaft`, `Generator`.
5. Reporting/plotting: live iteration tracking, result tables, performance map plots.
6. Multi-fluid / species tracking networks.
7. Transient simulation.
