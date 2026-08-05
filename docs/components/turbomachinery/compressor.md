# Compressor

<img src="../../_static/diagrams/compressor.svg" alt="Compressor diagram" class="component-diagram">

Two variants, same isentropic-compression physics: [`Compressor`](#compressor-map-driven)
is driven by a real characteristic map, [`SimpleCompressor`](#simplecompressor-fixed-ratio)
is an analytic fixed-ratio stand-in for when no map is available.

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters common to both:** `PR`
(`SimpleCompressor`) or a map (`Compressor`), `eta_s` (`SimpleCompressor`
only — `Compressor` reads it off the map), `gamma` (optional, else derived
from the fluid), `heat_path` (optional)

$$
T_\text{out,s} = T_\text{in}\,PR^{\frac{\gamma-1}{\gamma}}
\qquad
\Delta h_\text{actual} = \frac{h(P_\text{out}, T_\text{out,s}) - h_\text{in}}{\eta_s}
$$

$$
P_\text{out} - PR\cdot P_\text{in} = 0
\qquad
h_\text{out} - \left(h_\text{in} + \Delta h_\text{actual}\right) + \frac{Q_\text{loss}}{\dot m} = 0
$$

Inefficiency shows up as *extra* enthalpy rise (extra heating of the fluid)
— the standard adiabatic-irreversible-compression model. `heat_path` (a
`Convection`/`Conduction`/`Radiation` from
[Thermal network](../thermal-network/index.md)) is optional and defaults to
fully adiabatic (`Q_loss = 0`). Both variants share this relation; they only
differ in where `PR`/`eta_s` come from.

## `Compressor` (map-driven)

Driven by a Flownex-style (`.cop`) characteristic map instead of a fixed
`PR`/`eta_s` — pressure ratio and isentropic efficiency are read off the
map's iso-speed curves at the compressor's shaft speed `N` and the corrected
mass flow implied by the current inlet state.

**Additional ports:** mechanical: `shaft` &nbsp;·&nbsp; signal: `power`
&nbsp;·&nbsp; **Parameters:** `map_path`, `N` (shaft speed `[rev/min]`;
leave `None` to solve for it), `factor_overrides` (optional map
calibration)

$$
A = \frac{N/60}{\sqrt{T_\text{in}}} \quad\text{(corrected speed)}
\qquad
B = \frac{\dot m\sqrt{T_\text{in}}}{P_\text{in}} \quad\text{(corrected flow)}
$$

$$
PR = \text{map.pressure\_ratio}(A, B) \qquad \eta_s = \text{map.efficiency}(A, B)
$$

then the shared isentropic relation above using this `PR` and `eta_s`.
**Free parameter:** leave `N=None` to solve for shaft speed instead of a
fixed `PR`/power — `N` lives on the `"shaft"` mechanical port either way
(fixed when given, a free Newton unknown when `None`), needing a matching
[`Setpoint`](../control/setpoint.md)/[`Controller`](../control/controller.md)
(`free_param="shaft"`) or [`Shaft`](../mechanical/shaft.md) connection
(`kind="mechanical"`) to close it. `factor_overrides` lets you calibrate
the map's own conversion factors against test data without editing the map
file. The `"power"` signal port publishes the same value as
`report_metrics()["power [W]"]`, for a `Generator`/`Shaft` to read via
`kind="signal"` — connecting the `"shaft"` mechanical port to a `Shaft`
wires this automatically, see [`Shaft`](../mechanical/shaft.md).

## `SimpleCompressor` (fixed ratio)

<img src="../../_static/diagrams/simple_compressor.svg" alt="SimpleCompressor diagram" class="component-diagram">

Analytic compressor: fixed pressure ratio and isentropic efficiency, no
performance map — no mechanical/signal ports, since there's no shaft speed
to speak of.

**Parameters:** `PR`, `eta_s` given directly instead of read off a map.

For the same physics driven by a real characteristic map instead, see
[`Compressor`](#compressor-map-driven) above.

---
Part of [Turbomachinery](index.md).
