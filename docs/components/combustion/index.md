# Combustion

Two models of the same physical role — adding fuel energy to the flow —
at different levels of fidelity.

```{toctree}
:maxdepth: 1

combustor
```

| Component | Role |
|---|---|
| [`Combustor`](combustor.md#combustor-chemical-equilibrium) | Cantera chemical-equilibrium combustion products, real adiabatic-flame-temperature physics. |
| [`SimpleCombustor`](combustor.md#simplecombustor-fixed-lhv) | Fixed lower-heating-value (LHV) heat release, no combustion chemistry. |
