# SteamTurbine

<img src="../../_static/diagrams/steam_turbine.svg" alt="SteamTurbine diagram" style="max-width:100%">

The wet-steam-correct counterpart to [`SimpleTurbine`](simple-turbine.md):
expands a (possibly condensing) vapor via an entropy-based isentropic path,
since the ideal-gas `T_out = T_in·(1/PR)^((γ-1)/γ)` relation is wrong once
the expansion crosses into the two-phase dome. Requires a fluid exposing
`entropy_ph`/`enthalpy_ps` (`CoolPropFluid`).

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters:** exactly one of
`P_out` (absolute exhaust pressure) or `PR`, `eta_s` (default 0.85)

$$
s_\text{in} = s(P_\text{in}, h_\text{in})
\qquad
h_\text{out,s} = h(P_\text{out}, s_\text{in})
\qquad
\Delta h_\text{actual} = \eta_s\,(h_\text{in} - h_\text{out,s})
$$

`report_metrics()` exposes `x_out [-]` (exhaust quality) — a real
steam-turbine concern, since excessive wetness erodes the last-stage
blades.

---
Part of [Turbomachinery](index.md).
