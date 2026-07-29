# Turbine

<img src="../../_static/diagrams/turbine.svg" alt="Turbine diagram" style="max-width:100%">

Driven by a Flownex-style (`.tur`) characteristic map — same corrected
groups as [`Compressor`](compressor.md) (referenced to the turbine's own
inlet conditions). The map's pressure-ratio column is the turbine's
*expansion* ratio `P_in/P_out`, the opposite sense from the compressor map.

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters:** `map_path`, `gamma`
(optional), `N` (leave `None` to solve for it), `factor_overrides`
(optional), `heat_path` (optional)

Same corrected-parameter and residual shape as
[`SimpleTurbine`](simple-turbine.md), with `PR` and `eta_s` read off
`map.pressure_ratio(A, B)`/`map.efficiency(A, B)` instead of fixed. **Free
parameter:** leave `N=None` to solve for shaft speed from a target
elsewhere in the network (a target power, a [`Shaft`](../mechanical/shaft.md)
speed-tie, ...) instead of giving it directly.

---
Part of [Turbomachinery](index.md).
