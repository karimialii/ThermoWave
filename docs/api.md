# API reference

Each subpackage re-exports its public classes at the top level, so
`thermowave.core`, `thermowave.components`, `thermowave.fluids`, and
`thermowave.maps` are the entry points below — the per-module import paths
(e.g. `thermowave.components.compressor.Compressor`) still work and are
unchanged.

## `thermowave.core`

Network assembly, the Newton solver, and transient time-marching.

```{eval-rst}
.. automodule:: thermowave.core
```

## `thermowave.components`

Every network component (sources/sinks, flow elements, turbomachinery,
heat-transfer elements, controllers, ...), all sharing the `BaseComponent`
interface.

```{eval-rst}
.. automodule:: thermowave.components
```

## `thermowave.fluids`

Fluid property models — ideal gas, ideal-gas mixtures, and the optional
CoolProp/Cantera-backed real-fluid models.

```{eval-rst}
.. automodule:: thermowave.fluids
```

## `thermowave.maps`

Turbomachinery characteristic maps and torque-speed maps.

```{eval-rst}
.. automodule:: thermowave.maps
```
