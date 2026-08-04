# MultiPassHeatExchanger

<img src="../../_static/diagrams/multi_pass_heat_exchanger.svg" alt="MultiPassHeatExchanger diagram" style="max-width:100%">

`MultiPassHeatExchanger` is a thin, always-UA/NTU-mode subclass of
[`HeatExchanger`](heat-exchanger.md) — its constructor
(`MultiPassHeatExchanger(name, UA, PR_hot, PR_cold, n_passes=1,
arrangement=..., correlation=None)`) and behavior are unchanged, but the
full model (effectiveness-NTU relations, `n_passes`, `arrangement`,
`shell_and_tube`'s F-correction, `custom` correlations) is now documented
once, on the [`HeatExchanger`](heat-exchanger.md) page — see it for the
full math.

**Ports:** `hot_in`, `hot_out`, `cold_in`, `cold_out` &nbsp;·&nbsp;
**Parameters:** `UA`, `PR_hot`, `PR_cold`, `n_passes` (default 1),
`arrangement` (`"counterflow"` / `"parallel"` / `"crossflow"` /
`"shell_and_tube"` / `"custom"`), `correlation` (required iff
`arrangement="custom"`)

`HeatExchanger(UA=...)` is equivalent — use `MultiPassHeatExchanger`
directly when you want "this exchanger's rating comes from geometry, not a
given effectiveness" explicit at the call site (or an `isinstance()`
check); either way, no math is duplicated between the two classes.

---
Part of [Heat exchangers & phase change](index.md).
