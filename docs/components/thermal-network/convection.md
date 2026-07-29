# Convection

<img src="../../_static/diagrams/heat_path.svg" alt="Heat path diagram" style="max-width:100%">

A passive heat path between two `TemperatureSource` endpoints (`a`, `b`) —
no flow ports, no algebraic residuals of its own. Other components'
`state_derivative()`/`residuals()` pull `Q(state)` from it directly. `Q` is
positive when `a` is hotter than `b` (heat flowing `a → b`) — the same sign
convention `Pipe`'s own `heat_loss` uses.

**Ports:** none &nbsp;·&nbsp; **Parameters:** `a`, `b`, `h`
(heat-transfer coefficient), `A` (area)

$$Q = h\,A\,(T_a - T_b)$$

Covers both free and forced convection — they differ only in how `h` is
physically obtained (natural buoyancy-driven flow vs. an imposed flow over
a surface), not in this formula.

For conduction through a solid instead, see [`Conduction`](conduction.md);
for radiative exchange, see [`Radiation`](radiation.md).

---
Part of [Thermal network](index.md).
