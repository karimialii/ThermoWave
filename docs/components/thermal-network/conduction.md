# Conduction

<img src="../../_static/diagrams/heat_path.svg" alt="Heat path diagram" style="max-width:100%">

A passive heat path between two `TemperatureSource` endpoints (`a`, `b`) —
no flow ports, no algebraic residuals of its own. Same zero-port,
zero-residual, `a → b` sign convention as [`Convection`](convection.md).

**Ports:** none &nbsp;·&nbsp; **Parameters:** `a`, `b`, `k` (conductivity),
`A` (cross-sectional area), `L` (length)

$$Q = \frac{kA}{L}\,(T_a - T_b) = UA\,(T_a - T_b)$$

Steady 1D conduction through a solid path — e.g. a shaft conducting heat
between two casings.

---
Part of [Thermal network](index.md).
