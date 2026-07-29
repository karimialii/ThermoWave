# SimpleEvaporator

<img src="../../_static/diagrams/simple_evaporator.svg" alt="SimpleEvaporator diagram" style="max-width:100%">

Single-stream evaporator/boiler/superheater: adds heat to hit a specified
outlet condition without modeling the heat-source stream explicitly (the
phase-change analogue of [`SimpleCombustor`](../combustion/simple-combustor.md)).

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters:** `PR` (default 1.0),
`outlet_quality` (default 1.0 = saturated vapor), `superheat` (K above
saturation), `duty` (optional — fix heat directly instead of an outlet
spec)

- **Outlet-spec mode** (default, `duty=None`): outlet pinned to
  `outlet_quality`/`superheat` (same targets as [`Evaporator`](evaporator.md));
  `Q = mdot * (h_out - h_in)` is a *reported* result.
- **Duty mode** (`duty` given, `[W]`): `h_out = h_in + duty/mdot`; the
  resulting outlet quality/superheat is reported instead.

For the same physics coupled to an explicit heat-source stream, see
[`Evaporator`](evaporator.md).

---
Part of [Heat exchangers & phase change](index.md).
