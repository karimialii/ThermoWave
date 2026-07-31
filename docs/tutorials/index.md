# Tutorials

Longer, narrative guides that build a complete model end to end, including
the mistakes that are easy to make along the way and why the fix is what it
is.

```{toctree}
:maxdepth: 1

building-a-gas-turbine-model
building-a-combustion-chamber-model
```

- **[Building a gas-turbine model](building-a-gas-turbine-model.md)** — a
  single-shaft recuperated microturbine (T100-class), from an empty
  `Network` to a converged steady state to a PID-controlled transient. Covers
  the seven-step build shape, the thermal (casing/shaft) network, and four
  traps that cost real debugging time: where mechanical loss actually
  belongs, why more fuel can *lower* turbine outlet temperature at fixed
  load, the no-equilibrium fuel cliff, and why sweeps must run ascending.
- **[Building a 1D combustion-chamber model](building-a-combustion-chamber-model.md)**
  — a T100-class combustor split into primary (root) and dilution air,
  discretized into a liner of `Pipe`/`ThermalMass`/`Convection`/`Radiation`/
  `Conduction` stations instead of one lumped `Combustor`. Covers calibrating
  an unknown split fraction against measured data, and three traps: a
  too-small annulus diameter that looks like a convergence failure, why
  `Junction` doesn't reconcile mismatched branch pressures, and why the
  liner runs far hotter than a real metal survives without film cooling.

Looking for shorter, single-purpose snippets instead? See
[Examples](../examples/index.md). For what each component actually computes,
see the [Components](../components/index.md) reference.
