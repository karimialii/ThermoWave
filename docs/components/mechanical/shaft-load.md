# ShaftLoad

<img src="../../_static/diagrams/shaft_load.svg" alt="ShaftLoad diagram" style="max-width:100%">

A constant-power mechanical load on a [`Shaft`](shaft.md) — the electrical
demand a generator (plus power electronics) actually places on the rotor.
Unlike [`SimpleGenerator`](simple-generator.md)/[`Generator`](generator.md)
(passive readers that only *report* power), `ShaftLoad` is a genuine torque
on the shaft: its demanded power enters the shaft power balance directly,
so the equilibrium speed becomes wherever turbine power covers compressor
power *plus* this demand.

**Ports:** none &nbsp;·&nbsp; **Parameters:** `power` (electrical demand,
`[W]`), `efficiency` (default 1.0)

**Residuals:** none — it contributes no speed unknown, only a power term
`power / efficiency` to whichever `Shaft` lists it (with sign `-1.0` for a
draw). `power` is a plain mutable attribute, so a
[`Schedule`](../control/schedule.md) can drive it through a dispatch
profile during `solve_transient()`.

---
Part of [Mechanical & electrical](index.md).
