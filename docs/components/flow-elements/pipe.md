# Pipe

<img src="../../_static/diagrams/pipe.svg" alt="Pipe diagram" class="component-diagram">

Darcy-Weisbach friction pressure drop, with optional heat loss.
Discretized into `n_elem` sub-elements (default 1), each contributing a
momentum and an energy residual from the pipe's own constant inlet `mdot`.

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters:** `L` (length), `D`
(diameter), `roughness` `[m]` + `mu` (dynamic viscosity `[Pa*s]`) to compute
the Darcy friction factor from Reynolds number, or `f` given directly to
skip that calculation, `n_elem`, `heat_loss` `[W]` (optional, fixed, total
across the whole pipe), `heat_path` (optional
`Convection`/`Conduction`/`Radiation`, live-computed instead of fixed — see
below), `z_in`/`z_out` `[m]` (elevation, both `0.0` by default — see
below)

By default (no `f` given), the friction factor is recomputed every residual
call from the current Reynolds number `Re = ρvD/μ`:

$$
f = \begin{cases}
64 / Re & Re < 2300 \text{ (laminar)} \\
\text{Haaland or Colebrook-White} & Re \geq 2300 \text{ (turbulent)}
\end{cases}
$$

using roughness `ε` as the relative roughness `ε/D` in whichever turbulent
correlation is selected via `thermowave.core.settings.friction_correlation`
(`"haaland"` by default, explicit; or `"colebrook"`, an implicit fixed-point
iteration seeded from the Haaland value) — a project-wide default, not a
per-`Pipe` setting. Passing `f` directly instead (e.g. `f=0.0` for a
frictionless duct, a common way to model a heater/combustor segment via
`heat_loss` alone) skips this calculation entirely and uses the fixed value,
same as before.

**Per element** (length `L/n_elem`, heat `heat_loss/n_elem`):

$$
\rho = \rho(P, h) \qquad v = \frac{\dot m}{\rho A} \qquad A = \frac{\pi D^2}{4}
$$

$$
\Delta p_\text{friction} = f\,\frac{L_\text{elem}}{D}\,\frac{\rho v^2}{2}
\qquad\qquad
\Delta p_\text{gravity} = \rho\, g\, \Delta z_\text{elem}
$$

$$
P_\text{in} - P_\text{out} - \Delta p_\text{friction} - \Delta p_\text{gravity} = 0
\qquad\qquad
h_\text{in} - h_\text{out} - \frac{q_\text{elem}}{\dot m} = 0
$$

where `g = 9.80665 m/s^2` and `Δz_elem = (z_out - z_in) / n_elem`, using
that same element's own (inlet-based) `ρ` — no separate elevation model,
just one more term in the momentum residual already there. Rising elevation
(`z_out > z_in`) lowers `P_out` relative to `P_in` the way a real fluid
column's weight does, on top of friction; falling elevation raises it.
Leaving `z_in`/`z_out` at their `0.0` default makes this term exactly zero,
identical to every `Pipe` built before this parameter existed.

This is enough, with no other new machinery, to model **thermosiphon /
natural-circulation** loops: chain a heated riser (`z_out > z_in`) and a
cooled downcomer (`z_in > z_out`) into a closed loop with a free mass flow
(e.g. `Source(mdot=None)`) and no pump — the loop's existing momentum
residuals settle at whatever `mdot` balances the two legs' friction against
the buoyancy created by their density difference, the same way the solver
already picks any other free unknown. See
[A closed thermosiphon (natural-circulation) loop](../../examples/index.md#a-closed-thermosiphon-natural-circulation-loop)
for a full worked example.

plus one mass residual tying the outlet `mdot` back to the inlet's
(`q_elem` above is `heat_loss` and `heat_path`'s live `Q(state)` summed
together, then split evenly across `n_elem` — same accounting `Combustor`,
`Turbine`, and `Compressor` already use for their own `heat_path`).

Use `heat_path` (via `Network.add_heat_path()`) instead of the fixed
`heat_loss` when the loss should depend on live state — a wall temperature
that isn't known up front, or a radiative term that scales with `T^4`. A
common pattern: chain several single-element `Pipe`s to represent stations
along a duct, each with its own `Convection`/`Radiation` to its own
`ThermalMass` wall ring, and `Conduction` between neighboring rings — see
[Building a 1D combustion-chamber model](../../tutorials/building-a-combustion-chamber-model.md)
for a worked example.

---
Part of [Flow elements](index.md).
