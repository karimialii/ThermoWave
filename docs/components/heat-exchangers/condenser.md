# Condenser

<img src="../../_static/diagrams/condenser.svg" alt="Condenser diagram" style="max-width:100%">

Two-stream phase-change equipment: a condensing working fluid (`wf`)
coupled to an explicit single-phase coolant stream. Duty comes from the
`wf` outlet spec, **never** an effectiveness/UA calculation (`cp` is
effectively infinite during phase change).

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

For the same physics without modeling the coolant stream explicitly, see
[`SimpleCondenser`](simple-condenser.md); for the boiling counterpart, see
[`Evaporator`](evaporator.md).

---
Part of [Heat exchangers & phase change](index.md).
