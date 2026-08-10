# Drum

<img src="../../_static/diagrams/drum.svg" alt="Drum diagram" class="component-diagram">

A steam drum: a constant-volume, two-phase reservoir at saturation, with
genuine liquid-level dynamics — a boiler's steam/water separator.

**Ports:** `feed_in`, `steam_out`, `water_out` (plus `riser_in` if
`has_riser=True`, the default) &nbsp;·&nbsp; **Parameters:** `V`, `P0`,
`fluid`, `level0` (default 0.5), `has_riser`, `heat_loss` (optional),
`P_target`, `water_out_mdot`, `level_target` (all optional — see below)

**Differential state:** `(P, h)`, driven by mass and energy conservation on
the drum's contents (mass `m = ρ(P,h)·V`; the same finite-difference
2×2 linear solve as [`Tank`](tank.md), generalized to multiple inlets):

$$
\frac{dm}{dt} = \sum \dot m_\text{in} - \sum \dot m_\text{out}
\qquad
\frac{d(mu)}{dt} = \sum \dot m_\text{in} h_\text{in} - \sum \dot m_\text{out} h_\text{out} - \dot Q_\text{loss}
$$

Each outlet leaves at its own saturation enthalpy — `steam_out` at $h_g(P)$,
`water_out` at $h_f(P)$:

$$
P_\text{steam} - P_\text{drum} = 0 \qquad h_\text{steam} - h_g(P_\text{drum}) = 0
$$

$$
P_\text{water} - P_\text{drum} = 0 \qquad h_\text{water} - h_f(P_\text{drum}) = 0
$$

All mdots are left free by default, closed by whatever the drum connects
to. A drum's liquid level has **no steady-state restoring force** on its
own — it's a pure integrator, exactly like a real drum needs level control
— so a plain `Network.solve()` (`dt=None`) is singular in the drum's own
`h` unless `level_target` closes it (below); `Network.solve_transient()`
always works regardless, since `dt` drives `(P, h)` forward from `P0`/`h0`
directly.

**Closing the remaining degrees of freedom** (optional, add only the ones a
given topology actually needs — `Network.check_wiring()` names any that are
still missing):

- **`P_target`** `[Pa]`: adds a residual pinning the drum's own pressure
  directly, $P_\text{drum} - P_\text{target} = 0$ — a real boiler drum is
  essentially isobaric with its feed, so this needs no `Setpoint`/
  `Controller` detour to say so.
- **`water_out_mdot`** `[kg/s]`: adds a residual pinning `water_out`'s own
  mass flow directly (with `has_riser=True`, this is the riser/downcomer
  recirculation rate) — the same role `Source`'s own fixed `mdot` plays,
  for the same reason: nothing else conserves mass into a specific numeric
  value on its own.
- **`level_target`** `[-]`, e.g. `0.5` for a half-full drum: closes the
  steady-state singularity in `h` described above. Without it, `Network.
  solve()` doesn't just converge slowly — the Jacobian column for the
  drum's own `h` is exactly zero at the solution (`state_derivative()`'s
  own mass/energy balance reduces to `-h·(mdot_in - mdot_out)` in `h`,
  vanishing identically once mass balances), which shows up as the
  residual norm dropping to within ~1% of converging, then jumping back up
  by 2-3 orders of magnitude and oscillating indefinitely, rather than a
  clean singular-matrix error. Leave at the default `None` only for
  `Network.solve_transient()`, where `dt` already drives `h` forward and no
  steady closure is wanted.

---
Part of [Heat exchangers & phase change](index.md).
