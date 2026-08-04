# Turbine

<img src="../../_static/diagrams/turbine.svg" alt="Turbine diagram" style="max-width:100%">

Driven by a Flownex-style (`.tur`) characteristic map — same corrected
groups as [`Compressor`](compressor.md) (referenced to the turbine's own
inlet conditions). The map's pressure-ratio column is the turbine's
*expansion* ratio `P_in/P_out`, the opposite sense from the compressor map.

**Ports:** `in`, `out` &nbsp;·&nbsp; mechanical: `shaft` &nbsp;·&nbsp;
signal: `power` &nbsp;·&nbsp; **Parameters:** `map_path`, `gamma`
(optional), `N` (leave `None` to solve for it), `factor_overrides`
(optional), `heat_path` (optional)

Same corrected-parameter and residual shape as
[`SimpleTurbine`](simple-turbine.md), with `PR` and `eta_s` read off
`map.pressure_ratio(A, B)`/`map.efficiency(A, B)` instead of fixed. **Free
parameter:** leave `N=None` to solve for shaft speed from a target
elsewhere in the network (a target power via
[`Setpoint`](../control/setpoint.md) with `free_param="shaft"`, or a
[`Shaft`](../mechanical/shaft.md) speed-tie via
`network.connect(shaft, "m_i", turbine, "shaft", kind="mechanical")`)
instead of giving it directly. The `"power"` signal port publishes the same
value as `report_metrics()["power [W]"]`, for a `Generator`/`Shaft` to read
via `kind="signal"`.

---
Part of [Turbomachinery](index.md).
