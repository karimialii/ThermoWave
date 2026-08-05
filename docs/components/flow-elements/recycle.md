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

---
Part of [Flow elements](index.md).
