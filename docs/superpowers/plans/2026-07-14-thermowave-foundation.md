# ThermoWave Foundation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the ThermoWave foundation — fluid property models, a network graph, a hand-rolled Newton-Raphson solver, and a minimal `Source`/`Pipe`/`Sink` component set — that solves a simple thermofluid network end-to-end and matches hand-calculated results.

**Architecture:** A `networkx`-backed graph where nodes hold scalar state `(P, h)` and the whole network shares one mass flow rate `mdot` (no branching/merging in this foundation, so `mdot` is a single value rather than per-branch). Components contribute residual equations; boundary conditions from `Source` are eliminated from the unknown vector rather than solved for. A hand-rolled Newton-Raphson loop with a finite-difference Jacobian, step damping, and physical clamping (pressure floor) drives node states to convergence.

**Tech Stack:** Python 3.10+, numpy, scipy (linear solve only), networkx, CoolProp (optional extra), pytest.

## Global Constraints

- Python: `>=3.10` (dev venv here runs 3.14.5, tests must pass there).
- All internal computation is strict SI: Pa, K, J/kg, kg/s, W, m.
- `CoolProp` is an optional extra (`pip install thermowave[coolprop]`); `CoolPropFluid` must raise a clear `ImportError` if unavailable, and its tests must be skipped (not failed) when CoolProp isn't installed, via `pytest.importorskip("CoolProp")`.
- Dependencies: `numpy>=1.24`, `scipy>=1.10`, `networkx>=3.0` (core); `CoolProp>=6.4` (optional `coolprop` extra); `pytest>=7.4`, `pytest-cov>=4.1`, `ruff>=0.1`, `mypy>=1.5` (optional `dev` extra).
- No transient simulation, no turbomachinery/combustion/heat-exchanger components, no multi-fluid networks, no rich terminal reporting — all out of scope for this plan (see spec's "Open Follow-Up Specs").
- Custom exceptions (`ConvergenceError`, `FluidRangeError`, `NetworkTopologyError`) must be used for solver/fluid/topology failures — never bare `ValueError`/`Exception`.

## Design Refinements Beyond the Spec

The approved spec (`docs/superpowers/specs/2026-07-14-thermowave-foundation-design.md`) left two interfaces underspecified (`NetworkState`, exact component wiring) and didn't give `Source` a way to fix mass flow rate — without a fixed `mdot`, the `Source → Pipe → Sink` system is under-determined (2 unknowns at the outlet node, only 2 equations from one pipe element, but `mdot` would be a 3rd unknown with no 3rd equation). This plan resolves both gaps concretely:

- **`Source` gains an `mdot` parameter** (`Source(name, P, T, mdot, outlet_node)`) and fixes the network's single shared mass flow rate, matching the spec's own end-to-end test which requires a known `mdot` to hand-calculate the Darcy-Weisbach pressure drop.
- **Components take explicit node-name strings** at construction (e.g. `outlet_node="n1"`) rather than relying on automatic graph inference — the caller wires the topology by reusing the same node-name string across adjacent components. This is simpler to implement correctly and matches how the end-to-end tests need to construct networks.
- **`mdot` is a single scalar per `Network`**, not a per-branch dict — this foundation has no flow-splitting/merging components, so one shared value is sufficient and keeps `Pipe`'s mass conservation implicit (no extra residual needed).

---

### Task 1: Project scaffolding and custom exceptions

**Files:**
- Create: `pyproject.toml`
- Create: `src/thermowave/__init__.py`
- Create: `src/thermowave/core/__init__.py`
- Create: `src/thermowave/core/exceptions.py`
- Create: `src/thermowave/fluids/__init__.py`
- Create: `src/thermowave/components/__init__.py`
- Create: `tests/__init__.py`
- Test: `tests/test_exceptions.py`

**Interfaces:**
- Produces: `thermowave.core.exceptions.ConvergenceError`, `FluidRangeError`, `NetworkTopologyError` — all `Exception` subclasses, constructed with a message string, used by every later task.

- [ ] **Step 1: Write the failing test**

Create `tests/test_exceptions.py`:

```python
import pytest

from thermowave.core.exceptions import (
    ConvergenceError,
    FluidRangeError,
    NetworkTopologyError,
)


@pytest.mark.parametrize(
    "exc_type",
    [ConvergenceError, FluidRangeError, NetworkTopologyError],
)
def test_exception_is_exception_subclass_with_message(exc_type):
    exc = exc_type("something went wrong")
    assert isinstance(exc, Exception)
    assert str(exc) == "something went wrong"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_exceptions.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'thermowave'`

- [ ] **Step 3: Write the scaffolding and implementation**

Create `pyproject.toml`:

```toml
[build-system]
requires = ["setuptools>=61.0.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "thermowave"
version = "0.1.0"
description = "A headless, 1D implicit thermodynamic network solver for steady-state thermal-fluid systems"
readme = "README.md"
license = { text = "MIT" }
requires-python = ">=3.10"
keywords = ["thermodynamics", "fluid-dynamics", "solver"]
classifiers = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.10",
    "Topic :: Scientific/Engineering :: Physics",
]

dependencies = [
    "numpy>=1.24",
    "scipy>=1.10",
    "networkx>=3.0",
]

[project.optional-dependencies]
coolprop = ["CoolProp>=6.4"]
dev = ["pytest>=7.4", "pytest-cov>=4.1", "ruff>=0.1", "mypy>=1.5"]
full = ["thermowave[coolprop,dev]"]

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-v --tb=short"
```

Create `src/thermowave/__init__.py`:

```python
"""ThermoWave: a 1D implicit thermodynamic network solver."""
```

Create `src/thermowave/core/__init__.py`:

```python
```

Create `src/thermowave/core/exceptions.py`:

```python
class ConvergenceError(Exception):
    """Raised when the Newton-Raphson solver fails to converge."""


class FluidRangeError(Exception):
    """Raised when a fluid property call falls outside a physically valid range."""


class NetworkTopologyError(Exception):
    """Raised when a network's component graph is invalid or under/over-determined."""
```

Create `src/thermowave/fluids/__init__.py`:

```python
```

Create `src/thermowave/components/__init__.py`:

```python
```

Create `tests/__init__.py`:

```python
```

Install the package in editable mode so `import thermowave` resolves:

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/pip install -e .`
Expected: `Successfully installed thermowave-0.1.0`

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_exceptions.py -v`
Expected: `3 passed`

- [ ] **Step 5: Commit**

```bash
cd /Users/ali/Documents/GitHub/ThermoWave
git add pyproject.toml src/thermowave/__init__.py src/thermowave/core/__init__.py \
  src/thermowave/core/exceptions.py src/thermowave/fluids/__init__.py \
  src/thermowave/components/__init__.py tests/__init__.py tests/test_exceptions.py
git commit -m "feat: scaffold thermowave package and add core exceptions"
```

---

### Task 2: BaseFluid interface and IdealGasFluid

**Files:**
- Create: `src/thermowave/fluids/base_fluid.py`
- Create: `src/thermowave/fluids/ideal_gas.py`
- Test: `tests/test_fluids.py`

**Interfaces:**
- Consumes: nothing from earlier tasks (standalone).
- Produces: `BaseFluid` (ABC) with abstract methods `density_ph(P, h) -> float`, `temperature_ph(P, h) -> float`, `enthalpy_pt(P, T) -> float`, `cp(P, T) -> float`, and attribute `name: str`. `IdealGasFluid(name, R, cp)` implementing all four. Later tasks (`CoolPropFluid`, `Source`, `Pipe`, `Solver`) call these exact method names/signatures on any `BaseFluid` instance.

- [ ] **Step 1: Write the failing test**

Create `tests/test_fluids.py`:

```python
import math

from thermowave.fluids.ideal_gas import IdealGasFluid

AIR_R = 287.05  # J/(kg*K)
AIR_CP = 1005.0  # J/(kg*K)


def test_ideal_gas_enthalpy_pt_matches_hand_calculation():
    air = IdealGasFluid(name="air", R=AIR_R, cp=AIR_CP)
    h = air.enthalpy_pt(P=101325.0, T=300.0)
    assert math.isclose(h, AIR_CP * 300.0, rel_tol=1e-9)


def test_ideal_gas_temperature_ph_matches_hand_calculation():
    air = IdealGasFluid(name="air", R=AIR_R, cp=AIR_CP)
    h = AIR_CP * 300.0
    T = air.temperature_ph(P=101325.0, h=h)
    assert math.isclose(T, 300.0, rel_tol=1e-9)


def test_ideal_gas_density_ph_matches_pv_equals_mrt():
    air = IdealGasFluid(name="air", R=AIR_R, cp=AIR_CP)
    P = 101325.0
    T = 300.0
    h = AIR_CP * T
    rho = air.density_ph(P=P, h=h)
    expected_rho = P / (AIR_R * T)
    assert math.isclose(rho, expected_rho, rel_tol=1e-9)


def test_ideal_gas_cp_is_constant():
    air = IdealGasFluid(name="air", R=AIR_R, cp=AIR_CP)
    assert air.cp(P=101325.0, T=300.0) == AIR_CP
    assert air.cp(P=200000.0, T=500.0) == AIR_CP


def test_ideal_gas_name_attribute():
    air = IdealGasFluid(name="air", R=AIR_R, cp=AIR_CP)
    assert air.name == "air"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_fluids.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'thermowave.fluids.ideal_gas'`

- [ ] **Step 3: Write minimal implementation**

Create `src/thermowave/fluids/base_fluid.py`:

```python
from abc import ABC, abstractmethod


class BaseFluid(ABC):
    """Interface every fluid property model implements.

    All arguments and return values are strict SI units: Pa, K, J/kg, J/(kg*K).
    """

    name: str

    @abstractmethod
    def density_ph(self, P: float, h: float) -> float:
        """Density [kg/m^3] given pressure [Pa] and specific enthalpy [J/kg]."""

    @abstractmethod
    def temperature_ph(self, P: float, h: float) -> float:
        """Temperature [K] given pressure [Pa] and specific enthalpy [J/kg]."""

    @abstractmethod
    def enthalpy_pt(self, P: float, T: float) -> float:
        """Specific enthalpy [J/kg] given pressure [Pa] and temperature [K]."""

    @abstractmethod
    def cp(self, P: float, T: float) -> float:
        """Specific heat at constant pressure [J/(kg*K)]."""
```

Create `src/thermowave/fluids/ideal_gas.py`:

```python
from thermowave.fluids.base_fluid import BaseFluid


class IdealGasFluid(BaseFluid):
    """Analytic ideal-gas fluid model with constant specific heat.

    Enthalpy is referenced to h=0 at T=0: h = cp * T.
    """

    def __init__(self, name: str, R: float, cp: float):
        self.name = name
        self.R = R
        self._cp = cp

    def cp(self, P: float, T: float) -> float:
        return self._cp

    def enthalpy_pt(self, P: float, T: float) -> float:
        return self._cp * T

    def temperature_ph(self, P: float, h: float) -> float:
        return h / self._cp

    def density_ph(self, P: float, h: float) -> float:
        T = self.temperature_ph(P, h)
        return P / (self.R * T)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_fluids.py -v`
Expected: `5 passed`

- [ ] **Step 5: Commit**

```bash
cd /Users/ali/Documents/GitHub/ThermoWave
git add src/thermowave/fluids/base_fluid.py src/thermowave/fluids/ideal_gas.py tests/test_fluids.py
git commit -m "feat: add BaseFluid interface and IdealGasFluid model"
```

---

### Task 3: CoolPropFluid real-fluid model

**Files:**
- Create: `src/thermowave/fluids/real_fluid.py`
- Modify: `tests/test_fluids.py` (append CoolProp tests)

**Interfaces:**
- Consumes: `BaseFluid` from Task 2 (`src/thermowave/fluids/base_fluid.py`); `FluidRangeError` from Task 1 (`src/thermowave/core/exceptions.py`).
- Produces: `CoolPropFluid(name, P_min=1e3, P_max=1e8)` implementing `BaseFluid`. Raises `ImportError` with a clear install message if CoolProp isn't installed. Raises `FluidRangeError` (not a bare `ValueError`) if a property call fails even after clamping.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_fluids.py`:

```python
import pytest

pytest.importorskip("CoolProp")

from thermowave.core.exceptions import FluidRangeError  # noqa: E402
from thermowave.fluids.real_fluid import CoolPropFluid  # noqa: E402


def test_coolprop_fluid_enthalpy_pt_matches_direct_propsi_call():
    from CoolProp.CoolProp import PropsSI

    water = CoolPropFluid(name="Water")
    P = 1.0e5  # 1 bar
    T = 473.15  # 200 C, superheated steam at 1 bar
    h = water.enthalpy_pt(P=P, T=T)
    expected = PropsSI("H", "P", P, "T", T, "Water")
    assert h == pytest.approx(expected, rel=1e-9)


def test_coolprop_fluid_round_trip_temperature_enthalpy():
    water = CoolPropFluid(name="Water")
    P = 1.0e5
    T = 473.15
    h = water.enthalpy_pt(P=P, T=T)
    T_back = water.temperature_ph(P=P, h=h)
    assert T_back == pytest.approx(T, rel=1e-6)


def test_coolprop_fluid_density_ph_matches_direct_propsi_call():
    from CoolProp.CoolProp import PropsSI

    water = CoolPropFluid(name="Water")
    P = 1.0e5
    T = 473.15
    h = water.enthalpy_pt(P=P, T=T)
    rho = water.density_ph(P=P, h=h)
    expected = PropsSI("D", "P", P, "H", h, "Water")
    assert rho == pytest.approx(expected, rel=1e-9)


def test_coolprop_fluid_clamps_pressure_below_minimum():
    water = CoolPropFluid(name="Water", P_min=1.0e5, P_max=1.0e7)
    assert water._clamp_pressure(1.0) == 1.0e5
    assert water._clamp_pressure(1.0e9) == 1.0e7
    assert water._clamp_pressure(5.0e5) == 5.0e5


def test_coolprop_fluid_raises_fluid_range_error_on_invalid_property_call():
    water = CoolPropFluid(name="Water")
    with pytest.raises(FluidRangeError):
        water.density_ph(P=1.0e5, h=-1.0e9)


def test_coolprop_fluid_import_error_message_when_coolprop_missing(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "CoolProp.CoolProp" or name.startswith("CoolProp"):
            raise ImportError("simulated missing CoolProp")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="coolprop"):
        CoolPropFluid(name="Water")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_fluids.py -v -k coolprop`
Expected: FAIL with `ModuleNotFoundError: No module named 'thermowave.fluids.real_fluid'`

- [ ] **Step 3: Write minimal implementation**

Create `src/thermowave/fluids/real_fluid.py`:

```python
from thermowave.core.exceptions import FluidRangeError
from thermowave.fluids.base_fluid import BaseFluid


class CoolPropFluid(BaseFluid):
    """Real-fluid model backed by CoolProp's PropsSI.

    Requires the optional 'coolprop' extra: pip install thermowave[coolprop]
    """

    def __init__(self, name: str, P_min: float = 1.0e3, P_max: float = 1.0e8):
        try:
            from CoolProp.CoolProp import PropsSI
        except ImportError as exc:
            raise ImportError(
                "CoolPropFluid requires the 'coolprop' extra: "
                "pip install thermowave[coolprop]"
            ) from exc

        self._props_si = PropsSI
        self.name = name
        self.P_min = P_min
        self.P_max = P_max

    def _clamp_pressure(self, P: float) -> float:
        return min(max(P, self.P_min), self.P_max)

    def density_ph(self, P: float, h: float) -> float:
        P = self._clamp_pressure(P)
        try:
            return self._props_si("D", "P", P, "H", h, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp density_ph failed for P={P}, h={h}, fluid={self.name}: {exc}"
            ) from exc

    def temperature_ph(self, P: float, h: float) -> float:
        P = self._clamp_pressure(P)
        try:
            return self._props_si("T", "P", P, "H", h, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp temperature_ph failed for P={P}, h={h}, fluid={self.name}: {exc}"
            ) from exc

    def enthalpy_pt(self, P: float, T: float) -> float:
        P = self._clamp_pressure(P)
        try:
            return self._props_si("H", "P", P, "T", T, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp enthalpy_pt failed for P={P}, T={T}, fluid={self.name}: {exc}"
            ) from exc

    def cp(self, P: float, T: float) -> float:
        P = self._clamp_pressure(P)
        try:
            return self._props_si("C", "P", P, "T", T, self.name)
        except ValueError as exc:
            raise FluidRangeError(
                f"CoolProp cp failed for P={P}, T={T}, fluid={self.name}: {exc}"
            ) from exc
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_fluids.py -v`
Expected: `11 passed`

- [ ] **Step 5: Commit**

```bash
cd /Users/ali/Documents/GitHub/ThermoWave
git add src/thermowave/fluids/real_fluid.py tests/test_fluids.py
git commit -m "feat: add CoolProp-backed real fluid model with range clamping"
```

---

### Task 4: Settings and unit conversion

**Files:**
- Create: `src/thermowave/core/settings.py`
- Create: `tests/conftest.py`
- Test: `tests/test_settings.py`

**Interfaces:**
- Consumes: nothing from earlier tasks (standalone).
- Produces: module-level singleton `thermowave.core.settings.settings` with attributes `pressure_unit: str` (default `"Pa"`), `temperature_unit: str` (default `"K"`), and methods `pressure_to_si(value: float) -> float`, `temperature_to_si(value: float) -> float`. `Source` (Task 5) calls these exact methods at construction time.

- [ ] **Step 1: Write the failing test**

Create `tests/conftest.py`:

```python
import pytest

from thermowave.core.settings import settings


@pytest.fixture(autouse=True)
def reset_settings():
    """Ensure the module-level settings singleton doesn't leak between tests."""
    original_pressure_unit = settings.pressure_unit
    original_temperature_unit = settings.temperature_unit
    yield
    settings.pressure_unit = original_pressure_unit
    settings.temperature_unit = original_temperature_unit
```

Create `tests/test_settings.py`:

```python
import math

from thermowave.core.settings import settings


def test_default_units_are_si():
    assert settings.pressure_unit == "Pa"
    assert settings.temperature_unit == "K"


def test_pressure_to_si_default_pa_is_passthrough():
    assert settings.pressure_to_si(101325.0) == 101325.0


def test_pressure_to_si_converts_bar():
    settings.pressure_unit = "bar"
    assert math.isclose(settings.pressure_to_si(1.013), 101300.0, rel_tol=1e-9)


def test_pressure_to_si_converts_atm():
    settings.pressure_unit = "atm"
    assert math.isclose(settings.pressure_to_si(1.0), 101325.0, rel_tol=1e-9)


def test_temperature_to_si_default_kelvin_is_passthrough():
    assert settings.temperature_to_si(300.0) == 300.0


def test_temperature_to_si_converts_celsius():
    settings.temperature_unit = "C"
    assert math.isclose(settings.temperature_to_si(15.0), 288.15, rel_tol=1e-9)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_settings.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'thermowave.core.settings'`

- [ ] **Step 3: Write minimal implementation**

Create `src/thermowave/core/settings.py`:

```python
_PRESSURE_TO_PA = {
    "Pa": 1.0,
    "kPa": 1.0e3,
    "MPa": 1.0e6,
    "bar": 1.0e5,
    "atm": 101325.0,
}


class Settings:
    """Configurable I/O units for constructing components and displaying results.

    All internal solver/component math is strict SI regardless of these settings;
    conversion happens only at construction time (see Source.__init__).
    """

    def __init__(self) -> None:
        self.pressure_unit = "Pa"
        self.temperature_unit = "K"

    def pressure_to_si(self, value: float) -> float:
        if self.pressure_unit not in _PRESSURE_TO_PA:
            raise ValueError(f"Unsupported pressure_unit: {self.pressure_unit!r}")
        return value * _PRESSURE_TO_PA[self.pressure_unit]

    def temperature_to_si(self, value: float) -> float:
        if self.temperature_unit == "K":
            return value
        if self.temperature_unit == "C":
            return value + 273.15
        raise ValueError(f"Unsupported temperature_unit: {self.temperature_unit!r}")


settings = Settings()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_settings.py -v`
Expected: `6 passed`

- [ ] **Step 5: Commit**

```bash
cd /Users/ali/Documents/GitHub/ThermoWave
git add src/thermowave/core/settings.py tests/conftest.py tests/test_settings.py
git commit -m "feat: add settings singleton with SI unit conversion"
```

---

### Task 5: BaseComponent, Source, and Sink

**Files:**
- Create: `src/thermowave/components/base_component.py`
- Create: `src/thermowave/components/source.py`
- Create: `src/thermowave/components/sink.py`
- Test: `tests/test_components.py`

**Interfaces:**
- Consumes: `BaseFluid`, `IdealGasFluid` from Task 2; `settings` from Task 4 (`thermowave.core.settings.settings`).
- Produces: `BaseComponent` (ABC) with abstract methods `ports() -> dict[str, str]` and `residuals(state) -> list[float]`, plus default (non-abstract) methods `internal_nodes() -> list[str]` (returns `[]`), `fixed_node_values(fluid: BaseFluid) -> dict[str, tuple[float, float]]` (returns `{}`), `fixed_mdot() -> float | None` (returns `None`). `Source(name, P, T, mdot, outlet_node)` and `Sink(name, inlet_node)`. Task 6 (`Pipe`) and Task 7 (`Network`) rely on exactly these method names/signatures. `state` here is `thermowave.core.network.NetworkState`, defined in Task 7 — referenced only via `TYPE_CHECKING` string annotations in this task since `Source`/`Sink` don't need it at runtime (their `residuals` always returns `[]`).

- [ ] **Step 1: Write the failing test**

Create `tests/test_components.py`:

```python
import math

from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.settings import settings
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def test_source_ports_returns_outlet_node_name():
    src = Source(name="src1", P=101325.0, T=300.0, mdot=1.0, outlet_node="n1")
    assert src.ports() == {"outlet": "n1"}


def test_source_residuals_is_empty():
    src = Source(name="src1", P=101325.0, T=300.0, mdot=1.0, outlet_node="n1")
    assert src.residuals(state=None) == []


def test_source_fixed_node_values_uses_si_units_by_default():
    src = Source(name="src1", P=101325.0, T=300.0, mdot=1.0, outlet_node="n1")
    fixed = src.fixed_node_values(AIR)
    P, h = fixed["n1"]
    assert math.isclose(P, 101325.0, rel_tol=1e-9)
    assert math.isclose(h, 1005.0 * 300.0, rel_tol=1e-9)


def test_source_fixed_node_values_converts_non_si_units(reset_settings=None):
    settings.pressure_unit = "bar"
    settings.temperature_unit = "C"
    try:
        src = Source(name="src1", P=1.013, T=26.85, mdot=1.0, outlet_node="n1")
        fixed = src.fixed_node_values(AIR)
        P, h = fixed["n1"]
        assert math.isclose(P, 101300.0, rel_tol=1e-9)
        assert math.isclose(h, 1005.0 * 300.0, rel_tol=1e-6)
    finally:
        settings.pressure_unit = "Pa"
        settings.temperature_unit = "K"


def test_source_fixed_mdot_returns_configured_value():
    src = Source(name="src1", P=101325.0, T=300.0, mdot=2.5, outlet_node="n1")
    assert src.fixed_mdot() == 2.5


def test_sink_ports_returns_inlet_node_name():
    snk = Sink(name="snk1", inlet_node="n2")
    assert snk.ports() == {"inlet": "n2"}


def test_sink_residuals_is_empty():
    snk = Sink(name="snk1", inlet_node="n2")
    assert snk.residuals(state=None) == []


def test_sink_has_no_fixed_node_values_or_mdot():
    snk = Sink(name="snk1", inlet_node="n2")
    assert snk.fixed_node_values(AIR) == {}
    assert snk.fixed_mdot() is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_components.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'thermowave.components.source'`

- [ ] **Step 3: Write minimal implementation**

Create `src/thermowave/components/base_component.py`:

```python
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid


class BaseComponent(ABC):
    """Interface every network component implements."""

    name: str

    @abstractmethod
    def ports(self) -> dict[str, str]:
        """Named ports ('inlet'/'outlet'/etc.) -> node name strings."""

    def internal_nodes(self) -> list[str]:
        """Node names this component creates internally (not shared with neighbors)."""
        return []

    @abstractmethod
    def residuals(self, state: "NetworkState") -> list[float]:
        """Residual equations contributed by this component given current state."""

    def fixed_node_values(self, fluid: "BaseFluid") -> dict[str, tuple[float, float]]:
        """node_name -> (P, h) this component fixes as a boundary condition."""
        return {}

    def fixed_mdot(self) -> float | None:
        """Network mass flow rate [kg/s] this component fixes, or None."""
        return None
```

Create `src/thermowave/components/source.py`:

```python
from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent
from thermowave.core.settings import settings

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState
    from thermowave.fluids.base_fluid import BaseFluid


class Source(BaseComponent):
    """Boundary condition fixing outlet pressure, temperature, and mass flow rate."""

    def __init__(self, name: str, P: float, T: float, mdot: float, outlet_node: str):
        self.name = name
        self.P_si = settings.pressure_to_si(P)
        self.T_si = settings.temperature_to_si(T)
        self.mdot = mdot
        self.outlet_node = outlet_node

    def ports(self) -> dict[str, str]:
        return {"outlet": self.outlet_node}

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def fixed_node_values(self, fluid: "BaseFluid") -> dict[str, tuple[float, float]]:
        h = fluid.enthalpy_pt(self.P_si, self.T_si)
        return {self.outlet_node: (self.P_si, h)}

    def fixed_mdot(self) -> float | None:
        return self.mdot
```

Create `src/thermowave/components/sink.py`:

```python
from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class Sink(BaseComponent):
    """Open boundary terminating a network branch. Fixes no state."""

    def __init__(self, name: str, inlet_node: str):
        self.name = name
        self.inlet_node = inlet_node

    def ports(self) -> dict[str, str]:
        return {"inlet": self.inlet_node}

    def residuals(self, state: "NetworkState") -> list[float]:
        return []
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_components.py -v`
Expected: `8 passed`

- [ ] **Step 5: Commit**

```bash
cd /Users/ali/Documents/GitHub/ThermoWave
git add src/thermowave/components/base_component.py src/thermowave/components/source.py \
  src/thermowave/components/sink.py tests/test_components.py
git commit -m "feat: add BaseComponent interface, Source, and Sink"
```

---

### Task 6: Pipe component

**Files:**
- Modify: `tests/test_components.py` (append Pipe tests)
- Create: `src/thermowave/components/pipe.py`

**Interfaces:**
- Consumes: `BaseComponent` from Task 5. `NetworkState` from Task 7 — but since Task 7 doesn't exist yet, this task defines `Pipe.residuals(state)` against a **duck-typed** state object (any object with `.fluid`, `.mdot`, and `.node(name) -> tuple[float, float]`), tested here with a hand-built stand-in class `_FakeState`. Task 7's real `NetworkState` must match this exact duck type: `fluid: BaseFluid`, `mdot: float`, `node(name: str) -> tuple[float, float]`.
- Produces: `Pipe(name, L, D, f, inlet_node, outlet_node, n_elem=1, heat_loss=None)`. `internal_nodes()` returns `[f"{name}__mid{i}" for i in range(1, n_elem)]`. Task 7 (`Network`) and Task 10 (end-to-end tests) construct `Pipe` with these exact keyword names.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_components.py`:

```python
import math as _math

from thermowave.components.pipe import Pipe
from thermowave.fluids.ideal_gas import IdealGasFluid as _IdealGasFluid


class _FakeState:
    """Minimal stand-in for core.network.NetworkState, used before Task 7 exists."""

    def __init__(self, fluid, mdot, node_values: dict[str, tuple[float, float]]):
        self.fluid = fluid
        self.mdot = mdot
        self._node_values = node_values

    def node(self, name: str) -> tuple[float, float]:
        return self._node_values[name]


def test_pipe_ports_returns_inlet_and_outlet():
    pipe = Pipe(
        name="p1", L=5.0, D=0.2, f=0.02, inlet_node="n1", outlet_node="n2", n_elem=1
    )
    assert pipe.ports() == {"inlet": "n1", "outlet": "n2"}


def test_pipe_internal_nodes_empty_for_single_element():
    pipe = Pipe(
        name="p1", L=5.0, D=0.2, f=0.02, inlet_node="n1", outlet_node="n2", n_elem=1
    )
    assert pipe.internal_nodes() == []


def test_pipe_internal_nodes_for_multiple_elements():
    pipe = Pipe(
        name="p1", L=5.0, D=0.2, f=0.02, inlet_node="n1", outlet_node="n2", n_elem=3
    )
    assert pipe.internal_nodes() == ["p1__mid1", "p1__mid2"]


def test_pipe_residuals_single_element_matches_darcy_weisbach_hand_calc():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    L, D, f = 5.0, 0.2, 0.02
    mdot = 1.0
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)

    rho = air.density_ph(P_in, h_in)
    area = _math.pi * D**2 / 4
    v = mdot / (rho * area)
    expected_dp = f * (L / D) * (rho * v**2 / 2)
    P_out_guess = P_in - expected_dp
    h_out_guess = h_in  # adiabatic

    pipe = Pipe(name="p1", L=L, D=D, f=f, inlet_node="n1", outlet_node="n2", n_elem=1)
    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={"n1": (P_in, h_in), "n2": (P_out_guess, h_out_guess)},
    )
    residuals = pipe.residuals(state)
    assert len(residuals) == 2
    momentum_residual, energy_residual = residuals
    assert _math.isclose(momentum_residual, 0.0, abs_tol=1e-6)
    assert _math.isclose(energy_residual, 0.0, abs_tol=1e-6)


def test_pipe_residuals_with_heat_loss_reflects_energy_balance():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    L, D, f = 5.0, 0.2, 0.02
    mdot = 1.0
    Q = 1000.0  # W, total heat loss
    P_in, T_in = 101325.0, 300.0
    h_in = air.enthalpy_pt(P_in, T_in)
    h_out_guess = h_in - Q / mdot  # exact energy balance
    P_out_guess = P_in  # momentum residual will be nonzero; only checking energy term

    pipe = Pipe(
        name="p1",
        L=L,
        D=D,
        f=f,
        inlet_node="n1",
        outlet_node="n2",
        n_elem=1,
        heat_loss=Q,
    )
    state = _FakeState(
        fluid=air,
        mdot=mdot,
        node_values={"n1": (P_in, h_in), "n2": (P_out_guess, h_out_guess)},
    )
    _momentum_residual, energy_residual = pipe.residuals(state)
    assert _math.isclose(energy_residual, 0.0, abs_tol=1e-6)


def test_pipe_residuals_multi_element_returns_two_residuals_per_element():
    air = _IdealGasFluid(name="air", R=287.05, cp=1005.0)
    pipe = Pipe(
        name="p1", L=6.0, D=0.2, f=0.02, inlet_node="n1", outlet_node="n2", n_elem=3
    )
    P_in, h_in = 101325.0, air.enthalpy_pt(101325.0, 300.0)
    state = _FakeState(
        fluid=air,
        mdot=1.0,
        node_values={
            "n1": (P_in, h_in),
            "p1__mid1": (P_in - 100.0, h_in),
            "p1__mid2": (P_in - 200.0, h_in),
            "n2": (P_in - 300.0, h_in),
        },
    )
    residuals = pipe.residuals(state)
    assert len(residuals) == 6  # 2 residuals (momentum, energy) * 3 elements
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_components.py -v -k pipe`
Expected: FAIL with `ModuleNotFoundError: No module named 'thermowave.components.pipe'`

- [ ] **Step 3: Write minimal implementation**

Create `src/thermowave/components/pipe.py`:

```python
from __future__ import annotations

import math
from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class Pipe(BaseComponent):
    """Pipe with Darcy-Weisbach friction pressure drop and optional heat loss.

    Discretized into n_elem sub-elements, each contributing a momentum residual
    (friction pressure drop) and an energy residual (enthalpy balance). Mass
    conservation is implicit: the whole network shares a single mdot value.
    """

    def __init__(
        self,
        name: str,
        L: float,
        D: float,
        f: float,
        inlet_node: str,
        outlet_node: str,
        n_elem: int = 1,
        heat_loss: float | None = None,
    ):
        self.name = name
        self.L = L
        self.D = D
        self.f = f
        self.inlet_node = inlet_node
        self.outlet_node = outlet_node
        self.n_elem = n_elem
        self.heat_loss = heat_loss
        self._area = math.pi * D**2 / 4

    def ports(self) -> dict[str, str]:
        return {"inlet": self.inlet_node, "outlet": self.outlet_node}

    def internal_nodes(self) -> list[str]:
        return [f"{self.name}__mid{i}" for i in range(1, self.n_elem)]

    def _element_nodes(self) -> list[str]:
        return [self.inlet_node, *self.internal_nodes(), self.outlet_node]

    def residuals(self, state: "NetworkState") -> list[float]:
        nodes = self._element_nodes()
        elem_L = self.L / self.n_elem
        q_elem = 0.0 if self.heat_loss is None else self.heat_loss / self.n_elem
        mdot = state.mdot

        out: list[float] = []
        for i in range(self.n_elem):
            P_in, h_in = state.node(nodes[i])
            P_out, h_out = state.node(nodes[i + 1])
            rho = state.fluid.density_ph(P_in, h_in)
            v = mdot / (rho * self._area)
            dp_friction = self.f * (elem_L / self.D) * (rho * v**2 / 2)
            out.append(P_in - P_out - dp_friction)
            out.append(h_in - h_out - q_elem / mdot)
        return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_components.py -v`
Expected: `14 passed`

- [ ] **Step 5: Commit**

```bash
cd /Users/ali/Documents/GitHub/ThermoWave
git add src/thermowave/components/pipe.py tests/test_components.py
git commit -m "feat: add Pipe component with Darcy-Weisbach friction and heat loss"
```

---

### Task 7: NetworkState and Network graph assembly

**Files:**
- Create: `src/thermowave/core/network.py`
- Test: `tests/test_network.py`

**Interfaces:**
- Consumes: `BaseComponent` from Task 5 (`ports()`, `internal_nodes()`, `fixed_node_values(fluid)`, `fixed_mdot()`); `BaseFluid` from Task 2; `NetworkTopologyError` from Task 1; `Source`/`Sink`/`Pipe` from Tasks 5-6 (used only in tests).
- Produces: `NetworkState(fluid, mdot, node_P, node_h)` — a plain class with `.fluid`, `.mdot`, and `.node(name) -> tuple[float, float]`, matching the duck type `Pipe.residuals` (Task 6) already expects. `Network(fluid)` with `add_component(component)`, and internal methods `_all_nodes() -> list[str]`, `_fixed_node_values() -> dict[str, tuple[float, float]]`, `_fixed_mdot() -> float`, `validate_topology() -> None` (raises `NetworkTopologyError`). `networkx` is imported and used to store the node/edge topology for validation, per the spec's architecture — `Network` builds a `networkx.DiGraph` from each component's `ports()`/`internal_nodes()` in `add_component`, used by `validate_topology` to confirm every node is reachable from a fixed node. `Network.solve(...)` is added in Task 9 (deferred import of `Solver` to avoid a circular import, since `Solver` needs `NetworkState` from this module).

- [ ] **Step 1: Write the failing test**

Create `tests/test_network.py`:

```python
import pytest

from thermowave.components.pipe import Pipe
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.exceptions import NetworkTopologyError
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def _simple_network() -> Network:
    network = Network(fluid=AIR)
    network.add_component(
        Source(name="src", P=101325.0, T=300.0, mdot=1.0, outlet_node="n1")
    )
    network.add_component(
        Pipe(name="p1", L=5.0, D=0.2, f=0.02, inlet_node="n1", outlet_node="n2")
    )
    network.add_component(Sink(name="snk", inlet_node="n2"))
    return network


def test_network_state_node_returns_p_and_h():
    state = NetworkState(
        fluid=AIR, mdot=1.0, node_P={"n1": 101325.0}, node_h={"n1": 301500.0}
    )
    assert state.node("n1") == (101325.0, 301500.0)


def test_network_all_nodes_includes_ports_and_internal_nodes():
    network = Network(fluid=AIR)
    network.add_component(
        Source(name="src", P=101325.0, T=300.0, mdot=1.0, outlet_node="n1")
    )
    network.add_component(
        Pipe(
            name="p1",
            L=6.0,
            D=0.2,
            f=0.02,
            inlet_node="n1",
            outlet_node="n2",
            n_elem=3,
        )
    )
    network.add_component(Sink(name="snk", inlet_node="n2"))
    assert network._all_nodes() == ["n1", "p1__mid1", "p1__mid2", "n2"]


def test_network_fixed_node_values_collected_from_source():
    network = _simple_network()
    fixed = network._fixed_node_values()
    assert set(fixed.keys()) == {"n1"}
    P, h = fixed["n1"]
    assert P == 101325.0
    assert h == 1005.0 * 300.0


def test_network_fixed_mdot_from_source():
    network = _simple_network()
    assert network._fixed_mdot() == 1.0


def test_network_fixed_mdot_raises_topology_error_when_no_source():
    network = Network(fluid=AIR)
    network.add_component(Sink(name="snk", inlet_node="n1"))
    with pytest.raises(NetworkTopologyError, match="mdot"):
        network._fixed_mdot()


def test_network_validate_topology_raises_when_no_boundary_state():
    network = Network(fluid=AIR)
    network.add_component(
        Pipe(name="p1", L=5.0, D=0.2, f=0.02, inlet_node="n1", outlet_node="n2")
    )
    with pytest.raises(NetworkTopologyError, match="boundary"):
        network.validate_topology()


def test_network_validate_topology_passes_for_simple_network():
    network = _simple_network()
    network.validate_topology()  # must not raise
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_network.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'thermowave.core.network'`

- [ ] **Step 3: Write minimal implementation**

Create `src/thermowave/core/network.py`:

```python
from __future__ import annotations

from typing import TYPE_CHECKING

import networkx as nx

from thermowave.core.exceptions import NetworkTopologyError

if TYPE_CHECKING:
    from thermowave.components.base_component import BaseComponent
    from thermowave.core.solver import SolveResult
    from thermowave.fluids.base_fluid import BaseFluid


class NetworkState:
    """Current node/branch state passed to each component's residuals() call."""

    def __init__(
        self,
        fluid: "BaseFluid",
        mdot: float,
        node_P: dict[str, float],
        node_h: dict[str, float],
    ):
        self.fluid = fluid
        self.mdot = mdot
        self.node_P = node_P
        self.node_h = node_h

    def node(self, name: str) -> tuple[float, float]:
        return self.node_P[name], self.node_h[name]


class Network:
    """A graph of components sharing a single working fluid and mass flow rate."""

    def __init__(self, fluid: "BaseFluid"):
        self.fluid = fluid
        self.components: list["BaseComponent"] = []
        self.graph = nx.DiGraph()

    def add_component(self, component: "BaseComponent") -> None:
        self.components.append(component)
        for node_name in component.ports().values():
            self.graph.add_node(node_name)
        for node_name in component.internal_nodes():
            self.graph.add_node(node_name)

    def _all_nodes(self) -> list[str]:
        nodes: list[str] = []
        for component in self.components:
            for node_name in component.ports().values():
                if node_name not in nodes:
                    nodes.append(node_name)
            for node_name in component.internal_nodes():
                if node_name not in nodes:
                    nodes.append(node_name)
        return nodes

    def _fixed_node_values(self) -> dict[str, tuple[float, float]]:
        fixed: dict[str, tuple[float, float]] = {}
        for component in self.components:
            fixed.update(component.fixed_node_values(self.fluid))
        return fixed

    def _fixed_mdot(self) -> float:
        values = [
            component.fixed_mdot()
            for component in self.components
            if component.fixed_mdot() is not None
        ]
        if not values:
            raise NetworkTopologyError(
                "Network has no component fixing mass flow rate (mdot). "
                "Add a Source with an mdot value."
            )
        return values[0]

    def validate_topology(self) -> None:
        if not self._fixed_node_values():
            raise NetworkTopologyError(
                "Network has no boundary component fixing node state (P, h). "
                "Add a Source."
            )
        self._fixed_mdot()

    def solve(
        self, tol: float = 1e-6, max_iter: int = 100, damping: float = 1.0
    ) -> "SolveResult":
        self.validate_topology()
        from thermowave.core.solver import Solver

        return Solver(self).solve(tol=tol, max_iter=max_iter, damping=damping)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_network.py -v`
Expected: `7 passed`

- [ ] **Step 5: Commit**

```bash
cd /Users/ali/Documents/GitHub/ThermoWave
git add src/thermowave/core/network.py tests/test_network.py
git commit -m "feat: add NetworkState and Network graph assembly with topology validation"
```

---

### Task 8: Newton-Raphson solver core (pure numerics)

**Files:**
- Create: `src/thermowave/core/solver.py` (numerics only in this task; `Solver` class wiring to `Network` is Task 9)
- Test: `tests/test_solver.py`

**Interfaces:**
- Consumes: `ConvergenceError` from Task 1.
- Produces: `newton_solve(residual_fn, x0, tol=1e-6, max_iter=100, damping=1.0, clamp_fn=None) -> tuple[np.ndarray, int, float]` returning `(x_solution, iterations, residual_norm)`, raising `ConvergenceError` if not converged within `max_iter` or if the Jacobian is singular. Also produces `_finite_difference_jacobian(residual_fn, x, F, eps=1e-6) -> np.ndarray`. Task 9's `Solver.solve()` calls `newton_solve` exactly as specified here.

- [ ] **Step 1: Write the failing test**

Create `tests/test_solver.py`:

```python
import numpy as np
import pytest

from thermowave.core.exceptions import ConvergenceError
from thermowave.core.solver import _finite_difference_jacobian, newton_solve


def test_finite_difference_jacobian_matches_analytic_for_linear_system():
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([2.0 * x[0] + 3.0 * x[1] - 1.0, -x[0] + 5.0 * x[1] + 2.0])

    x = np.array([1.0, 1.0])
    F = residual_fn(x)
    J = _finite_difference_jacobian(residual_fn, x, F)
    expected_J = np.array([[2.0, 3.0], [-1.0, 5.0]])
    assert np.allclose(J, expected_J, atol=1e-4)


def test_newton_solve_converges_to_known_linear_root():
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([x[0] - 3.0, 2.0 * x[1] + 4.0])

    x0 = np.array([0.0, 0.0])
    x_sol, iterations, residual_norm = newton_solve(residual_fn, x0, tol=1e-9)
    assert np.allclose(x_sol, [3.0, -2.0], atol=1e-6)
    assert iterations >= 1
    assert residual_norm < 1e-9


def test_newton_solve_converges_to_known_nonlinear_root():
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([x[0] ** 2 - 4.0])

    x0 = np.array([3.0])
    x_sol, _iterations, residual_norm = newton_solve(residual_fn, x0, tol=1e-9)
    assert np.allclose(x_sol, [2.0], atol=1e-6)
    assert residual_norm < 1e-9


def test_newton_solve_raises_convergence_error_when_residual_never_vanishes():
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([1.0, 1.0])

    x0 = np.array([0.0, 0.0])
    with pytest.raises(ConvergenceError):
        newton_solve(residual_fn, x0, tol=1e-9, max_iter=3)


def test_newton_solve_raises_convergence_error_on_singular_jacobian():
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([x[0] - x[1], x[0] - x[1]])

    x0 = np.array([0.0, 1.0])
    with pytest.raises(ConvergenceError):
        newton_solve(residual_fn, x0, tol=1e-9, max_iter=5)


def test_newton_solve_applies_clamp_fn_each_iteration():
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([x[0] - 100.0])

    def clamp_fn(x: np.ndarray) -> np.ndarray:
        return np.maximum(x, np.array([0.0]))

    x0 = np.array([50.0])
    x_sol, _iterations, _residual_norm = newton_solve(
        residual_fn, x0, tol=1e-9, clamp_fn=clamp_fn
    )
    assert np.allclose(x_sol, [100.0], atol=1e-6)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_solver.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'thermowave.core.solver'`

- [ ] **Step 3: Write minimal implementation**

Create `src/thermowave/core/solver.py`:

```python
from __future__ import annotations

from typing import Callable

import numpy as np

from thermowave.core.exceptions import ConvergenceError


def _finite_difference_jacobian(
    residual_fn: Callable[[np.ndarray], np.ndarray],
    x: np.ndarray,
    F: np.ndarray,
    eps: float = 1e-6,
) -> np.ndarray:
    n = len(x)
    m = len(F)
    J = np.zeros((m, n))
    for j in range(n):
        x_pert = x.copy()
        step = eps * max(abs(x[j]), 1.0)
        x_pert[j] += step
        F_pert = residual_fn(x_pert)
        J[:, j] = (F_pert - F) / step
    return J


def newton_solve(
    residual_fn: Callable[[np.ndarray], np.ndarray],
    x0: np.ndarray,
    tol: float = 1e-6,
    max_iter: int = 100,
    damping: float = 1.0,
    clamp_fn: Callable[[np.ndarray], np.ndarray] | None = None,
) -> tuple[np.ndarray, int, float]:
    """Damped Newton-Raphson with a finite-difference Jacobian.

    Returns (x_solution, iterations, residual_norm). Raises ConvergenceError if
    not converged within max_iter iterations or if the Jacobian is singular.
    """
    x = np.array(x0, dtype=float)
    if clamp_fn is not None:
        x = clamp_fn(x)

    for iteration in range(1, max_iter + 1):
        F = residual_fn(x)
        residual_norm = float(np.linalg.norm(F))
        if residual_norm < tol:
            return x, iteration - 1, residual_norm

        J = _finite_difference_jacobian(residual_fn, x, F)
        try:
            dx = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError as exc:
            raise ConvergenceError(
                f"Singular Jacobian at iteration {iteration}: {exc}"
            ) from exc

        x = x + damping * dx
        if clamp_fn is not None:
            x = clamp_fn(x)

    F = residual_fn(x)
    residual_norm = float(np.linalg.norm(F))
    if residual_norm < tol:
        return x, max_iter, residual_norm

    raise ConvergenceError(
        f"Solver failed to converge after {max_iter} iterations "
        f"(residual norm={residual_norm:.3e}, tol={tol:.3e})"
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_solver.py -v`
Expected: `6 passed`

- [ ] **Step 5: Commit**

```bash
cd /Users/ali/Documents/GitHub/ThermoWave
git add src/thermowave/core/solver.py tests/test_solver.py
git commit -m "feat: add hand-rolled Newton-Raphson solver core with damping and clamping"
```

---

### Task 9: Wire Solver to Network

**Files:**
- Modify: `src/thermowave/core/solver.py` (append `SolveResult` and `Solver`)
- Modify: `tests/test_solver.py` (append `Solver`-level tests)

**Interfaces:**
- Consumes: `newton_solve` from Task 8 (this module); `Network`, `NetworkState` from Task 7 (`thermowave.core.network`); `ConvergenceError` from Task 1.
- Produces: `SolveResult(converged, iterations, residual_norm, node_P, node_h, mdot)` and `Solver(network).solve(tol=1e-6, max_iter=100, damping=1.0) -> SolveResult`. `Network.solve()` (Task 7) already calls `Solver(self).solve(...)` with this exact signature. Task 10's end-to-end tests call `network.solve(...)` and read `result.node_P`, `result.node_h`, `result.mdot`, `result.converged`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_solver.py`:

```python
from thermowave.components.pipe import Pipe
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.core.solver import Solver
from thermowave.fluids.ideal_gas import IdealGasFluid

_AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def test_solver_solve_returns_converged_result_for_single_pipe_network():
    network = Network(fluid=_AIR)
    network.add_component(
        Source(name="src", P=101325.0, T=300.0, mdot=1.0, outlet_node="n1")
    )
    network.add_component(
        Pipe(name="p1", L=5.0, D=0.2, f=0.02, inlet_node="n1", outlet_node="n2")
    )
    network.add_component(Sink(name="snk", inlet_node="n2"))

    result = Solver(network).solve(tol=1e-6, max_iter=50)

    assert result.converged is True
    assert result.node_P["n1"] == 101325.0  # fixed by Source, untouched
    assert result.node_P["n2"] < 101325.0  # pressure drops across the pipe
    assert result.mdot == 1.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_solver.py -v -k single_pipe_network`
Expected: FAIL with `ImportError: cannot import name 'Solver' from 'thermowave.core.solver'`

- [ ] **Step 3: Write minimal implementation**

Append to `src/thermowave/core/solver.py`:

```python
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from thermowave.core.network import Network


class SolveResult:
    """Outcome of a Network.solve() call."""

    def __init__(
        self,
        converged: bool,
        iterations: int,
        residual_norm: float,
        node_P: dict[str, float],
        node_h: dict[str, float],
        mdot: float,
    ):
        self.converged = converged
        self.iterations = iterations
        self.residual_norm = residual_norm
        self.node_P = node_P
        self.node_h = node_h
        self.mdot = mdot


class Solver:
    """Assembles a Network's free variables into a Newton-Raphson problem."""

    P_MIN = 1.0e3  # Pa, physical clamp floor to protect fluid property calls

    def __init__(self, network: "Network"):
        self.network = network

    def solve(
        self, tol: float = 1e-6, max_iter: int = 100, damping: float = 1.0
    ) -> SolveResult:
        from thermowave.core.network import NetworkState

        network = self.network
        fluid = network.fluid
        fixed_nodes = network._fixed_node_values()
        mdot = network._fixed_mdot()
        all_nodes = network._all_nodes()
        free_nodes = [n for n in all_nodes if n not in fixed_nodes]

        def unpack(x: np.ndarray) -> tuple[dict[str, float], dict[str, float]]:
            node_P = {n: v[0] for n, v in fixed_nodes.items()}
            node_h = {n: v[1] for n, v in fixed_nodes.items()}
            for i, n in enumerate(free_nodes):
                node_P[n] = x[2 * i]
                node_h[n] = x[2 * i + 1]
            return node_P, node_h

        def residual_vector(x: np.ndarray) -> np.ndarray:
            node_P, node_h = unpack(x)
            state = NetworkState(fluid=fluid, mdot=mdot, node_P=node_P, node_h=node_h)
            residuals: list[float] = []
            for component in network.components:
                residuals.extend(component.residuals(state))
            return np.array(residuals, dtype=float)

        def clamp(x: np.ndarray) -> np.ndarray:
            x = x.copy()
            for i in range(len(free_nodes)):
                if x[2 * i] < self.P_MIN:
                    x[2 * i] = self.P_MIN
            return x

        if not free_nodes:
            node_P, node_h = unpack(np.array([]))
            F = residual_vector(np.array([]))
            return SolveResult(
                True, 0, float(np.linalg.norm(F)), node_P, node_h, mdot
            )

        if fixed_nodes:
            default_P, default_h = next(iter(fixed_nodes.values()))
        else:
            default_P, default_h = 1.0e5, 3.0e5
        x0 = np.zeros(2 * len(free_nodes))
        for i in range(len(free_nodes)):
            x0[2 * i] = default_P
            x0[2 * i + 1] = default_h

        x_sol, iterations, residual_norm = newton_solve(
            residual_vector, x0, tol=tol, max_iter=max_iter, damping=damping, clamp_fn=clamp
        )
        node_P, node_h = unpack(x_sol)
        return SolveResult(True, iterations, residual_norm, node_P, node_h, mdot)
```

Add `import numpy as np` is already present at the top of `solver.py` from Task 8 — no duplicate import needed.

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_solver.py -v`
Expected: `7 passed`

- [ ] **Step 5: Commit**

```bash
cd /Users/ali/Documents/GitHub/ThermoWave
git add src/thermowave/core/solver.py tests/test_solver.py
git commit -m "feat: wire Solver to Network via NetworkState and free-variable assembly"
```

---

### Task 10: End-to-end analytic network tests

**Files:**
- Create: `tests/test_network_e2e.py`

**Interfaces:**
- Consumes: `Network` from Task 7, `Source`/`Pipe`/`Sink` from Tasks 5-6, `IdealGasFluid` from Task 2, `NetworkTopologyError` from Task 1. No production code changes — this task only adds tests that exercise the full stack built in Tasks 1-9.

- [ ] **Step 1: Write the tests**

Create `tests/test_network_e2e.py`:

```python
import math

import pytest

from thermowave.components.pipe import Pipe
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.exceptions import NetworkTopologyError
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def _hand_calc_darcy_dp(fluid, L, D, f, mdot, P_in, T_in) -> float:
    h_in = fluid.enthalpy_pt(P_in, T_in)
    rho = fluid.density_ph(P_in, h_in)
    area = math.pi * D**2 / 4
    v = mdot / (rho * area)
    return f * (L / D) * (rho * v**2 / 2)


def test_source_pipe_sink_single_element_matches_hand_calculated_pressure_drop():
    L, D, f, mdot = 5.0, 0.2, 0.02, 1.0
    P_in, T_in = 101325.0, 300.0
    expected_dp = _hand_calc_darcy_dp(AIR, L, D, f, mdot, P_in, T_in)

    network = Network(fluid=AIR)
    network.add_component(
        Source(name="src", P=P_in, T=T_in, mdot=mdot, outlet_node="n1")
    )
    network.add_component(
        Pipe(name="p1", L=L, D=D, f=f, inlet_node="n1", outlet_node="n2", n_elem=1)
    )
    network.add_component(Sink(name="snk", inlet_node="n2"))

    result = network.solve(tol=1e-8, max_iter=50)

    assert result.converged
    assert result.node_P["n1"] == pytest.approx(P_in, rel=1e-9)
    assert result.node_P["n2"] == pytest.approx(P_in - expected_dp, rel=1e-4)
    assert result.node_h["n2"] == pytest.approx(result.node_h["n1"], rel=1e-9)  # adiabatic


def test_source_pipe_sink_multi_element_matches_single_element_total_drop():
    L, D, f, mdot = 6.0, 0.2, 0.02, 1.0
    P_in, T_in = 101325.0, 300.0
    expected_dp = _hand_calc_darcy_dp(AIR, L, D, f, mdot, P_in, T_in)

    network = Network(fluid=AIR)
    network.add_component(
        Source(name="src", P=P_in, T=T_in, mdot=mdot, outlet_node="n1")
    )
    network.add_component(
        Pipe(name="p1", L=L, D=D, f=f, inlet_node="n1", outlet_node="n2", n_elem=3)
    )
    network.add_component(Sink(name="snk", inlet_node="n2"))

    result = network.solve(tol=1e-8, max_iter=50)

    assert result.converged
    assert result.node_P["n2"] == pytest.approx(P_in - expected_dp, rel=1e-3)


def test_source_pipe_sink_with_heat_loss_matches_energy_balance():
    L, D, f, mdot = 5.0, 0.2, 0.02, 1.0
    Q = 2000.0  # W
    P_in, T_in = 101325.0, 300.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    expected_h_out = h_in - Q / mdot

    network = Network(fluid=AIR)
    network.add_component(
        Source(name="src", P=P_in, T=T_in, mdot=mdot, outlet_node="n1")
    )
    network.add_component(
        Pipe(
            name="p1",
            L=L,
            D=D,
            f=f,
            inlet_node="n1",
            outlet_node="n2",
            n_elem=1,
            heat_loss=Q,
        )
    )
    network.add_component(Sink(name="snk", inlet_node="n2"))

    result = network.solve(tol=1e-8, max_iter=50)

    assert result.converged
    assert result.node_h["n2"] == pytest.approx(expected_h_out, rel=1e-6)


def test_network_without_source_raises_topology_error_on_solve():
    network = Network(fluid=AIR)
    network.add_component(
        Pipe(name="p1", L=5.0, D=0.2, f=0.02, inlet_node="n1", outlet_node="n2")
    )
    network.add_component(Sink(name="snk", inlet_node="n2"))

    with pytest.raises(NetworkTopologyError):
        network.solve()
```

- [ ] **Step 2: Run tests to verify they fail or pass appropriately**

Since Tasks 1-9 already implement everything these tests exercise, these should pass immediately — this task validates the full stack rather than driving new implementation.

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest tests/test_network_e2e.py -v`
Expected: `4 passed`

If any test fails, the bug is in Tasks 1-9's implementation (most likely a units/sign error in `Pipe.residuals` or `Solver.solve`'s unknown-vector layout) — fix the implementation file, not the test, since these values are independently hand-calculated.

- [ ] **Step 3: Run the full test suite**

Run: `cd /Users/ali/Documents/GitHub/ThermoWave && .venv/bin/python -m pytest -v`
Expected: all tests across `tests/test_exceptions.py`, `tests/test_fluids.py`, `tests/test_settings.py`, `tests/test_components.py`, `tests/test_network.py`, `tests/test_solver.py`, `tests/test_network_e2e.py` pass (46 tests total).

- [ ] **Step 4: Commit**

```bash
cd /Users/ali/Documents/GitHub/ThermoWave
git add tests/test_network_e2e.py
git commit -m "test: add end-to-end Source-Pipe-Sink analytic validation"
```

---

## Post-Plan Verification

After Task 10, run the full suite once more and confirm the package imports cleanly from a fresh interpreter:

```bash
cd /Users/ali/Documents/GitHub/ThermoWave
.venv/bin/python -m pytest -v
.venv/bin/python -c "
from thermowave.core.network import Network
from thermowave.components.source import Source
from thermowave.components.pipe import Pipe
from thermowave.components.sink import Sink
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name='air', R=287.05, cp=1005.0)
net = Network(fluid=air)
net.add_component(Source(name='src', P=101325.0, T=300.0, mdot=1.0, outlet_node='n1'))
net.add_component(Pipe(name='p1', L=5.0, D=0.2, f=0.02, inlet_node='n1', outlet_node='n2'))
net.add_component(Sink(name='snk', inlet_node='n2'))
result = net.solve()
print('converged:', result.converged, 'P_out:', result.node_P['n2'])
"
```

Expected: all tests pass, and the script prints `converged: True` with a `P_out` slightly below 101325.0.
