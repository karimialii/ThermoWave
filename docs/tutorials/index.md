# Tutorials

Longer, narrative guides that build a complete model end to end, including
the mistakes that are easy to make along the way and why the fix is what it
is.

```{toctree}
:maxdepth: 1

building-a-gas-turbine-model
building-a-combustion-chamber-model
closing-a-recycle-composition-loop
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
  `Conduction` stations instead of one lumped `Combustor`. Covers deriving
  the split from an independent design rule (not fitting it to measured
  data), an engineering trade study that finds a genuine temperature-vs-CO
  conflict with no single satisfying split, why equilibrium NOx numbers
  aren't usable as compliance predictions, and three traps: a too-small
  annulus diameter that looks like a convergence failure, why `Junction`
  doesn't reconcile mismatched branch pressures, and why the liner runs far
  hotter than a real metal survives without film cooling.
- **[Closing a recycle/EGR composition loop](closing-a-recycle-composition-loop.md)**
  — why a `Recycle`-closed flow loop and a *composition* cycle through that
  same loop are two different problems, and why the second one needs
  `Recycle(fluid_guess=...)`. Covers the tear-stream (direct-substitution)
  iteration `Network.solve()` runs to converge it, a worked fresh-air/
  recirculated-stream `Junction` example with a checkable analytic answer,
  and two traps: wiring `Recycle` directly across a component with a forced
  `(P, h)` change (over-determined, never converges), and why composition
  doesn't propagate through a second consecutive `Combustor` pass.

Looking for shorter, single-purpose snippets instead? See
[Examples](../examples/index.md). For what each component actually computes,
see the [Components](../components/index.md) reference.
