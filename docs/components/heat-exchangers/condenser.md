# Condenser

<img src="../../_static/diagrams/condenser.svg" alt="Condenser diagram" class="component-diagram">

Two variants: [`Condenser`](#condenser-two-stream) couples the condensing
working fluid to an explicit coolant stream, and
[`SimpleCondenser`](#simplecondenser-single-stream) hits the same outlet
condition without modeling the coolant. Both are phase-change equipment:
duty always comes from the outlet spec, **never** an effectiveness/UA
calculation (`cp` is effectively infinite during phase change).

## `Condenser` (two-stream)

A condensing working fluid (`wf`) coupled to an explicit single-phase
coolant stream.

**Ports:** `wf_in`, `wf_out`, `cool_in`, `cool_out` &nbsp;·&nbsp;
**Parameters:** `PR_wf`, `PR_cool`, `outlet_quality` (default 0.0 =
saturated liquid), `subcool` (K below saturation, if > 0)

$$
h_\text{wf,out,target} = \begin{cases}
h(P_\text{wf,out}, q=\text{outlet\_quality}) & \text{subcool} = 0 \\
h\!\left(P_\text{wf,out},\, T_\text{sat}(P_\text{wf,out}) - \text{subcool}\right) & \text{subcool} > 0
\end{cases}
$$

$$
Q = \dot m_\text{wf}\cdot\left(h_\text{wf,in} - h_\text{wf,out,target}\right)
\qquad
h_\text{cool,out} = h_\text{cool,in} + \frac{Q}{\dot m_\text{cool}}
$$

**Pinch is a diagnostic, not a solved constraint**: `report_metrics()`
exposes `pinch [K]` `= T_sat(P_wf_out) - T_cool_out`. A negative pinch
means the coolant would have to end up hotter than the condensing fluid —
a thermodynamically infeasible spec, not a solver bug.

For the boiling counterpart, see [`Evaporator`](evaporator.md).

## `SimpleCondenser` (single-stream)

<img src="../../_static/diagrams/simple_condenser.svg" alt="SimpleCondenser diagram" class="component-diagram">

Rejects heat to hit a specified outlet condition without modeling the
coolant stream explicitly (the phase-change analogue of
[`SimpleCombustor`](../combustion/combustor.md)).

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters:** `PR` (default 1.0),
`outlet_quality` (default 0.0 = saturated liquid), `subcool` (K below
saturation), `duty` (optional — fix heat directly instead of an outlet
spec)

- **Outlet-spec mode** (default, `duty=None`): outlet pinned to
  `outlet_quality`/`subcool` (same targets as `Condenser` above);
  `Q = mdot * (h_out - h_in)` (negative) is a *reported* result.
- **Duty mode** (`duty` given, `[W]`, negative for heat rejection):
  `h_out = h_in + duty/mdot`; the resulting outlet quality/subcool is
  reported instead.

---
Part of [Heat exchangers & phase change](index.md).
