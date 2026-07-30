# Pipe: roughness-based friction factor

## Problem

`Pipe` currently requires a Darcy friction factor `f` given directly as a constant. Real pipe friction depends on pipe roughness and the actual flow's Reynolds number (which varies with mdot, fluid properties, and pipe diameter) — a fixed `f` doesn't reflect that physics.

## Goals

- Compute the Darcy friction factor from pipe roughness and Reynolds number by default, re-evaluated fresh every residual call (mdot, and hence Re, changes during the Newton solve).
- Cover both laminar (`f = 64/Re`) and turbulent (roughness-dependent) regimes.
- Support two turbulent correlations — Haaland (explicit) and Colebrook-White (implicit, fixed-point) — selectable via a new global default on the existing `thermowave.core.settings` singleton, defaulting to Haaland.
- Keep `f` as an optional direct override (same pattern as `Turbine.gamma`): when given, skip the roughness/Re calculation entirely. This keeps every existing `f=0.02`/`f=0.0` call site working unchanged, including the ~15 tests/docs that use `f=0.0` deliberately to model a frictionless heater duct.

## Non-goals

- Deriving viscosity from the fluid model. No `BaseFluid` implementation in this codebase exposes viscosity today; `Pipe` takes `mu` as a direct constant, the same trust level `gamma` already has on Turbine/Compressor.
- A per-Pipe correlation override — the choice is a single project-wide default on `settings`, not a per-Pipe kwarg (per explicit user choice).
- Smoothly blending laminar/turbulent across the transition Reynolds number — a hard switch at Re=2300 is acceptable for this codebase's finite-difference-Jacobian Newton solver.

## Design

### 1. `thermowave/core/settings.py`

Add `self.friction_correlation = "haaland"` to `Settings.__init__`, validated the same way `pressure_unit`/`temperature_unit` already are (raise `ValueError` on an unsupported value), against `{"haaland", "colebrook"}`.

### 2. `thermowave/components/pipe.py`

Constructor signature:
```python
def __init__(
    self,
    name: str,
    L: float,
    D: float,
    roughness: float | None = None,
    mu: float | None = None,
    f: float | None = None,
    n_elem: int = 1,
    heat_loss: float | None = None,
):
```
- If `f` is `None`, both `roughness` and `mu` must be given (`ValueError` otherwise).
- Store `roughness`, `mu`, `f` as given.

Friction factor, computed inside `residuals()` per element (roughness/D and mu are constant along the pipe; only rho/v change per element):
```python
def _friction_factor(self, rho: float, v: float) -> float:
    if self.f is not None:
        return self.f
    Re = max(rho * abs(v) * self.D / self.mu, 1e-9)
    if Re < 2300.0:
        return 64.0 / Re
    rel_roughness = self.roughness / self.D
    x_haaland = -1.8 * math.log10(
        (rel_roughness / 3.7) ** 1.11 + 6.9 / Re
    )
    if settings.friction_correlation == "haaland":
        return 1.0 / x_haaland**2
    # Colebrook-White: fixed-point iteration seeded from the Haaland value.
    x = x_haaland
    for _ in range(15):
        x = -2.0 * math.log10(rel_roughness / 3.7 + 2.51 / (Re * x))
    return 1.0 / x**2
```
`dp_friction = self._friction_factor(rho, v) * (elem_L / self.D) * (rho * v**2 / 2)` replaces the current `self.f * ...` line.

### 3. Docs

- `docs/components/flow-elements/pipe.md`: document `roughness`, `mu`, and the two-regime correlation; keep `f` documented as the direct-override escape hatch.
- The three canonical "single flow branch" examples (`README.md`, `docs/index.md`, `docs/examples/index.md`) switch their `Pipe(...)` call to `roughness=`/`mu=` to showcase the new physical path, since these are the first examples a new user reads.
- All other existing `Pipe(..., f=...)` call sites (tests, other docs) are left unchanged — they exercise the override path, which is intentionally unaffected.

## Testing

- Unit tests in `tests/test_pipe_friction.py`:
  - `f` given directly still bypasses the calculation (same dp as before).
  - Missing both `f` and `roughness`/`mu` raises `ValueError`.
  - Laminar branch: low Re gives `f = 64/Re` matching a hand-computed value.
  - Turbulent Haaland branch: matches a hand-computed Haaland value at a known Re/roughness.
  - Turbulent Colebrook branch: matches a reference Colebrook-White value (e.g. from a known Moody-chart point) within a tight tolerance, and differs from Haaland's value by a small amount (not identical, not wildly different).
  - `settings.friction_correlation` actually switches which branch a turbulent-regime pipe uses.
  - Invalid `settings.friction_correlation` value raises `ValueError`.
- Full existing suite must still pass unchanged (regression check on the override path).
