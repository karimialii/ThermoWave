# Component `set()` method — design

## Problem

Component constructors take many explicit kwargs (e.g. `Turbine(name, map_path, gamma=None, N=None, factor_overrides=None, heat_path=None)`). Setting several parameters inline in the constructor call makes it hard to scan a script later and tell which parameters were deliberately configured. Today the workaround is direct attribute assignment after construction (e.g. `turb.heat_path = path`), which already works (attributes are plain, unvalidated instance attributes — no dataclasses, no `**kwargs`) but has no typo protection and doesn't invalidate the owning `Network`'s caches.

## Goals

- A `component.set(**kwargs)` method usable any time before `Network.solve()`/`solve_transient()` — not restricted to "right after construction" — matching the existing pattern where physical parameters (`gamma`, `N`, `heat_path`, ...) are read live every residual/solver call.
- Reject unknown parameter names (typo protection).
- No per-component declaration of "settable" parameters — inherited once on `BaseComponent`, works for all 33 concrete components with zero changes to their files.
- Auto-invalidate the owning `Network`'s caches, removing the need to remember `network.invalidate_caches()` after a `set()` call.

## Non-goals

- Re-validating parameter semantics (e.g. Turbine's `gamma > 1` check) inside `set()`. Direct attribute assignment doesn't validate today either; `set()` matches that trust level.
- Preventing structurally unwise sets (e.g. overwriting `map`). No per-component allow-list is introduced.

## Design

### 1. `BaseComponent.set()` — `src/thermowave/components/base_component.py`

```python
def set(self, **kwargs) -> "BaseComponent":
    for key, value in kwargs.items():
        if key.startswith("_") or key == "name" or not hasattr(self, key):
            raise AttributeError(
                f"{type(self).__name__} {self.name!r} has no settable parameter {key!r}"
            )
        setattr(self, key, value)
    network = getattr(self, "_network", None)
    if network is not None:
        network.invalidate_caches()
    return self
```

- Validates against attributes that already exist on the instance post-`__init__` (`hasattr`) — this is what gives typo protection without a separate declared list.
- Excludes `name` (identity, not tunable) and any `_`-prefixed attribute (internal bookkeeping like `_inlet_node`).
- Returns `self` to allow chaining/fluent calls.

### 2. Back-reference in `Network.add_component()` — `src/thermowave/core/network.py:210`

```python
def add_component(self, component: "BaseComponent") -> None:
    self.components.append(component)
    component._network = self
    self._register_nodes(component)
```

This is the only change needed outside `base_component.py`. Every concrete component inherits `set()` for free.

### 3. Docs

Update docstrings that currently point at direct attribute assignment as the way to configure a component post-construction (e.g. Turbine's "just set `turb.heat_path = path` afterwards") to mention `turb.set(heat_path=path)` as the preferred spelling, keeping direct assignment mentioned as still working identically.

## Testing

- Unit test: `set()` on a fresh component updates the attribute and returns `self`.
- Unit test: `set()` with an unknown kwarg raises `AttributeError`.
- Unit test: `set()` with `name=...` or a `_`-prefixed key raises `AttributeError`.
- Unit test: after `network.add_component(turb)`, calling `turb.set(...)` bumps the network's topology version (i.e. `invalidate_caches()` was actually invoked), and calling `set()` on a component never added to a network doesn't error (no `_network` attribute).
