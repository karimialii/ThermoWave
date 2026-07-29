# Tank

<img src="../../_static/diagrams/tank.svg" alt="Tank diagram" style="max-width:100%">

A constant-volume plenum with real mass/energy storage — the genuinely
transient counterpart to [`Junction`](../flow-elements/junction.md)'s
zero-volume, quasi-steady mixing.

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters:** `V`, `P0`, `T0`,
`fluid`, `heat_loss` (optional)

**Differential state:** `(P, h)`, from the same mass/energy conservation as
[`Drum`](drum.md) (single inlet/outlet):

$$
\frac{dm}{dt} = \dot m_\text{in} - \dot m_\text{out}
\qquad
\frac{d(mu)}{dt} = \dot m_\text{in} h_\text{in} - \dot m_\text{out} h_\text{out} - \dot Q_\text{loss}
$$

solved as a 2×2 linear system for `dP/dt`, `dh/dt` using finite-differenced
`density_ph()` partials — works with any `BaseFluid`, not just calorically-
perfect gases. The outlet always carries the tank's own current well-mixed
state:

$$
P_\text{out} - P_\text{tank} = 0 \qquad h_\text{out} - h_\text{tank} = 0
$$

Neither `mdot` is constrained by the tank itself — real inflow and outflow
generally differ (that's the point of accumulation); each is closed by
whatever's upstream/downstream. `Network.solve()` solves for the `(P, h)` at
which both derivatives are exactly zero (a genuine steady-state plenum
condition); `Network.solve_transient()` integrates forward in time,
producing real filling/blow-down/surge behavior.

---
Part of [Heat exchangers & phase change](index.md).
