# Combustion

Two models of the same physical role — adding fuel energy to the flow —
at different levels of fidelity.

```{toctree}
:maxdepth: 1

simple-combustor
combustor
```

| Component | Role |
|---|---|
| [`SimpleCombustor`](simple-combustor.md) | Fixed lower-heating-value (LHV) heat release, no combustion chemistry. |
| [`Combustor`](combustor.md) | Cantera chemical-equilibrium combustion products, real adiabatic-flame-temperature physics. |
