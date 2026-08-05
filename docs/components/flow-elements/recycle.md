# Recycle

Mass-flow anchor for a genuinely closed loop — no `Source`/`Sink` needed.

**Ports:** `in`, `out`

**Fixes:** `mdot` at `out` only — never `(P, h)`. This is deliberately the
opposite half of [`Source`](source.md)'s job: `Source` always fixes `(P, h)`
and optionally `mdot`; `Recycle` always fixes `mdot` and never touches
`(P, h)`.

**Residuals:** 2 — an ideal, lossless pass-through:

$$P_\text{out} - P_\text{in} = 0 \qquad h_\text{out} - h_\text{in} = 0$$

Not the usual 3-residual flow-component pattern: there's deliberately no
mass-conservation residual, since `fixed_node_mdot()` already fixes this
component's own outlet `mdot` directly.

**Why this exists:** every flow component's own mass residual
(`mdot_out - mdot_in = 0`) already chains `mdot` all the way around a loop —
what's missing without a `Source` is an anchor for the shared *magnitude*,
not a relationship between the loop's mdot unknowns (conservation alone only
pins mdot up to an arbitrary constant). Likewise, `(P, h)` at every node in a
real closed cycle (ORC, Rankine, Brayton) are already fully determined by the
loop's own component physics — a [`SimpleCondenser`](../heat-exchangers/condenser.md)
computes its own outlet state from a subcool/quality target, a `Pump` raises
`P` by a given ratio, and so on — so nothing needs a hard boundary condition
either. `Recycle` supplies exactly the missing piece: drop it anywhere in a
loop (typically condenser outlet → pump inlet) to close the ring:

```python
net.connect(condenser, "out", recycle, "in")
net.connect(recycle, "out", pump, "in")
```

Use [`Source`](source.md)/[`Sink`](sink.md) for an open system with a real
inlet/exhaust boundary; use `Recycle` for a working fluid that recirculates.

## Composition cycles (EGR-style loops)

The flow (`P`/`h`/`mdot`) side above is handled by the ordinary Newton
solve regardless of whether the graph has a cycle in it. Fluid
**composition** is resolved by a separate pass that only ever propagates
forward from the network's fixed `(P, h)` boundaries — and a `Recycle`-closed
loop fixes none of those, so a node whose composition depends on itself
through the loop (a recirculated [`Combustor`](../combustion/combustor.md)
exhaust, a [`Junction`](junction.md) blend feeding back on itself) would
otherwise never resolve, silently falling back to the network's plain
default fluid.

`Recycle(fluid_guess=<a BaseFluid>)` seeds that pass with a starting
composition at its own outlet. `Network.solve()` then re-solves and
compares what actually arrives back at this component's inlet against the
guess (`recycle_fluid_delta()`), updates it (`update_fluid_guess()`), and
re-solves again — direct-substitution ("tear stream") iteration — until they
agree within `recycle_fluid_tol` (default `1e-4`), or raises
`ConvergenceError` after `recycle_fluid_max_iter` outer passes (default
`25`):

```python
recycle = Recycle(name="recycle", mdot=1.0, fluid_guess=some_fluid)
result = network.solve(recycle_fluid_tol=1e-6, recycle_fluid_max_iter=30)
```

Leave `fluid_guess` unset (the default) for a loop with no
composition-changing component in it — the ordinary Rankine/Brayton case
above needs none of this. See the
[recycle/EGR composition loop tutorial](../../tutorials/closing-a-recycle-composition-loop.md)
for a full worked example and its two traps.

---
Part of [Flow elements](index.md).
