# Turbine

<img src="../../_static/diagrams/turbine.svg" alt="Turbine diagram" class="component-diagram">

Three variants: [`Turbine`](#turbine-map-driven) is driven by a real
characteristic map, [`SimpleTurbine`](#simpleturbine-fixed-ratio) is an
analytic fixed-ratio stand-in for when no map is available, and
[`SteamTurbine`](#steamturbine-wet-steam) is the wet-steam-correct variant
for expansion through the two-phase dome.

**Ports common to all three:** `in`, `out`

## `Turbine` (map-driven)

Driven by a Flownex-style (`.tur`) characteristic map — same corrected
groups as [`Compressor`](compressor.md) (referenced to the turbine's own
inlet conditions). The map's pressure-ratio column is the turbine's
*expansion* ratio `P_in/P_out`, the opposite sense from the compressor map.

**Additional ports:** mechanical: `shaft` &nbsp;·&nbsp; signal: `power`
&nbsp;·&nbsp; **Parameters:** `map_path`, `gamma` (optional), `N` (leave
`None` to solve for it), `factor_overrides` (optional), `heat_path`
(optional)

Same corrected-parameter and residual shape as
[`SimpleTurbine`](#simpleturbine-fixed-ratio) below, with `PR` and `eta_s`
read off `map.pressure_ratio(A, B)`/`map.efficiency(A, B)` instead of
fixed. **Free parameter:** leave `N=None` to solve for shaft speed from a
target elsewhere in the network (a target power via
[`Setpoint`](../control/setpoint.md) with `free_param="shaft"`, or a
[`Shaft`](../mechanical/shaft.md) speed-tie via
`network.connect(shaft, "m_i", turbine, "shaft", kind="mechanical")`)
instead of giving it directly. The `"power"` signal port publishes the same
value as `report_metrics()["power [W]"]`, for a `Generator`/`Shaft` to read
via `kind="signal"` — connecting the `"shaft"` mechanical port to a `Shaft`
wires this automatically, see [`Shaft`](../mechanical/shaft.md).

## `SimpleTurbine` (fixed ratio)

<img src="../../_static/diagrams/simple_turbine.svg" alt="SimpleTurbine diagram" class="component-diagram">

Analytic turbine: fixed expansion ratio (`PR = P_in/P_out > 1`) and
isentropic efficiency, no performance map — no mechanical/signal ports,
since there's no shaft speed to speak of.

**Parameters:** `PR`, `eta_s`, `gamma` (optional), `heat_path` (optional)

$$
T_\text{out,s} = T_\text{in}\left(\frac{1}{PR}\right)^{\frac{\gamma-1}{\gamma}}
\qquad
\Delta h_\text{actual} = \eta_s\,\left(h_\text{in} - h(P_\text{out}, T_\text{out,s})\right)
$$

$$
P_\text{in} - PR\cdot P_\text{out} = 0
\qquad
h_\text{out} - \left(h_\text{in} - \Delta h_\text{actual}\right) + \frac{Q_\text{loss}}{\dot m} = 0
$$

Inefficiency means less work is extracted than the ideal expansion would
give.

## `SteamTurbine` (wet steam)

<img src="../../_static/diagrams/steam_turbine.svg" alt="SteamTurbine diagram" class="component-diagram">

The wet-steam-correct variant: expands a (possibly condensing) vapor via an
entropy-based isentropic path, since the ideal-gas
`T_out = T_in·(1/PR)^((γ-1)/γ)` relation `SimpleTurbine` uses above is
wrong once the expansion crosses into the two-phase dome. Requires a fluid
exposing `entropy_ph`/`enthalpy_ps` (`CoolPropFluid`).

**Parameters:** exactly one of `P_out` (absolute exhaust pressure) or `PR`,
`eta_s` (default 0.85)

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
