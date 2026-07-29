# ThermoWave

A headless, 1D implicit thermodynamic network solver for steady-state
thermal-fluid systems. You build a network out of components (sources, pipes,
valves, compressors, turbines, heat exchangers, combustors, shafts, sensors,
controllers, ...), wire their ports together, and call `solve()`. Internally
every component contributes unknowns and residual equations to one square
system, solved by damped Newton-Raphson with a finite-difference Jacobian.

See the project [README](https://github.com/karimialii/ThermoWave#readme)
for a quick-start example and installation instructions.

```{toctree}
:maxdepth: 2
:caption: Contents

building-a-gas-turbine-model
api
```
