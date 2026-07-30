# Thermal network

Solid thermal mass and the heat paths between them — casings, shafts, and
ambient, coupled to the flow network via a component's optional
`heat_path` attribute (see [Turbomachinery](../turbomachinery/index.md) and
[Combustion](../combustion/index.md)). None of these have flow ports; they
only participate in the differential-state/Newton bookkeeping, or are read
as a passive `Q(state)` by another component's own residual.

Every heat path's `a`/`b` endpoint is one of: a [`ThermalMass`](thermal-mass.md),
a fixed `float` (e.g. ambient temperature), or a `(component, port_name)`
tuple reading that component's live fluid temperature at that port's node.

```{toctree}
:maxdepth: 1

thermal-mass
convection
conduction
radiation
```

| Component | Role |
|---|---|
| [`ThermalMass`](thermal-mass.md) | A solid's own temperature as time-integrated state (a casing, a shaft, ...). |
| [`Convection`](convection.md) | `Q = h·A·(T_a - T_b)` — free or forced convection. |
| [`Conduction`](conduction.md) | `Q = (k·A/L)·(T_a - T_b)` — steady 1D conduction through a solid. |
| [`Radiation`](radiation.md) | `Q = ε·F·σ·A·(T_a⁴ - T_b⁴)` — surface-to-surface or surface-to-ambient radiative exchange. |

## Wiring a heat path into the flow network

[`Turbine`](../turbomachinery/turbine.md), [`Compressor`](../turbomachinery/compressor.md),
[`SimpleTurbine`](../turbomachinery/simple-turbine.md),
[`SimpleCompressor`](../turbomachinery/simple-compressor.md),
[`SimpleCombustor`](../combustion/simple-combustor.md), and
[`Combustor`](../combustion/combustor.md) each accept an optional
`heat_path` attribute (`None` by default — fully adiabatic). Since a path
needs `(component, "out")` as one of its own endpoints, it can only be
built *after* that component exists:

```python
from thermowave.components import Turbine, ThermalMass, Convection

turb = Turbine(name="turb", map_path="turbine.tur", N=65_000.0)
casing = ThermalMass(name="turb_casing", thermal_capacitance=200.0, T0=300.0)
conv = Convection(name="turb_conv", a=(turb, "out"), b=casing, h=50.0, A=0.3)
turb.set(heat_path=conv)      # or plain assignment, turb.heat_path = conv

network.add_heat_path(conv)   # registers the path and derives casing's sign
```

`Network.add_heat_path()` wires a path into both of its endpoints in one
call — deriving each endpoint's sign from whether it's `a` (loses `Q`) or
`b` (gains `Q`), and registering any `ThermalMass` endpoint so its
temperature doesn't go missing from the Newton system. A flow component can
carry a list of several heat paths at once (e.g. a combustor liner both
radiating to its casing and convecting to the annulus air), not just one.
