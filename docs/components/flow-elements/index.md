# Flow elements

Boundary conditions and passive flow-restriction elements — the pipes,
valves, and boundaries that connect the active machinery in
[Turbomachinery](../turbomachinery/index.md), [Combustion](../combustion/index.md),
and [Heat exchangers](../heat-exchangers/index.md) into a complete network.

```{toctree}
:maxdepth: 1

source
sink
pipe
valve
check-valve
nozzle
junction
```

| Component | Role |
|---|---|
| [`Source`](source.md) | Fixes outlet pressure, temperature, and (usually) mass flow — a network's inlet boundary. |
| [`Sink`](sink.md) | Terminates a branch, optionally pinning exit pressure. |
| [`Pipe`](pipe.md) | Darcy-Weisbach friction pressure drop, with optional heat loss. |
| [`Valve`](valve.md) | Flow-restriction pressure drop via a loss coefficient, isenthalpic. |
| [`CheckValve`](check-valve.md) | One-way variant of `Valve` — stiff resistance against reverse flow. |
| [`Nozzle`](nozzle.md) | Converts a pressure drop into velocity, with a choked-flow mass-flow cap. |
| [`Junction`](junction.md) | Splits or merges any number of flow branches. |
