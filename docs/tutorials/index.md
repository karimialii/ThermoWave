# Tutorials

Longer, narrative guides that build a complete model end to end, including
the mistakes that are easy to make along the way and why the fix is what it
is.

```{toctree}
:maxdepth: 1

building-a-gas-turbine-model
```

- **[Building a gas-turbine model](building-a-gas-turbine-model.md)** — a
  single-shaft recuperated microturbine (T100-class), from an empty
  `Network` to a converged steady state to a PID-controlled transient. Covers
  the seven-step build shape, the thermal (casing/shaft) network, and four
  traps that cost real debugging time: where mechanical loss actually
  belongs, why more fuel can *lower* turbine outlet temperature at fixed
  load, the no-equilibrium fuel cliff, and why sweeps must run ascending.

Looking for shorter, single-purpose snippets instead? See
[Examples](../examples/index.md). For what each component actually computes,
see the [Components](../components/index.md) reference.
