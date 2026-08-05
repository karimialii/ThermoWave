# Evaporator

<img src="../../_static/diagrams/evaporator.svg" alt="Evaporator diagram" class="component-diagram">

Two variants: [`Evaporator`](#evaporator-two-stream) couples the boiling
working fluid to an explicit heat-source stream, and
[`SimpleEvaporator`](#simpleevaporator-single-stream) hits the same outlet
condition without modeling the heat source. Both are phase-change
equipment: duty always comes from the outlet spec, **never** an
effectiveness/UA calculation (`cp` is effectively infinite during phase
change).

## `Evaporator` (two-stream)

A boiling working fluid (`wf`) coupled to an explicit single-phase
heat-source stream.

**Ports:** `wf_in`, `wf_out`, `src_in`, `src_out` &nbsp;·&nbsp;
**Parameters:** `PR_wf`, `PR_src`, `outlet_quality` (default 1.0 =
saturated vapor), `superheat` (K above saturation, if > 0)

$$
h_\text{wf,out,target} = \begin{cases}
h(P_\text{wf,out}, q=\text{outlet\_quality}) & \text{superheat} = 0 \\
h\!\left(P_\text{wf,out},\, T_\text{sat}(P_\text{wf,out}) + \text{superheat}\right) & \text{superheat} > 0
\end{cases}
$$

$$
Q = \dot m_\text{wf}\cdot\left(h_\text{wf,out,target} - h_\text{wf,in}\right)
\qquad
h_\text{src,out} = h_\text{src,in} - \frac{Q}{\dot m_\text{src}}
$$

**Pinch is a diagnostic, not a solved constraint**: `report_metrics()`
exposes `pinch [K]` `= T_src_out - T_sat(P_wf_out)`. A negative pinch means
the source would have to end up colder than the boiling fluid — a
thermodynamically infeasible spec, not a solver bug.

For the condensing counterpart, see [`Condenser`](condenser.md).

## `SimpleEvaporator` (single-stream)

<img src="../../_static/diagrams/simple_evaporator.svg" alt="SimpleEvaporator diagram" class="component-diagram">

Single-stream evaporator/boiler/superheater: adds heat to hit a specified
outlet condition without modeling the heat-source stream explicitly (the
phase-change analogue of [`SimpleCombustor`](../combustion/combustor.md)).

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters:** `PR` (default 1.0),
`outlet_quality` (default 1.0 = saturated vapor), `superheat` (K above
saturation), `duty` (optional — fix heat directly instead of an outlet
spec)

- **Outlet-spec mode** (default, `duty=None`): outlet pinned to
  `outlet_quality`/`superheat` (same targets as `Evaporator` above);
  `Q = mdot * (h_out - h_in)` is a *reported* result.
- **Duty mode** (`duty` given, `[W]`): `h_out = h_in + duty/mdot`; the
  resulting outlet quality/superheat is reported instead.

---
Part of [Heat exchangers & phase change](index.md).
