# Compressor

<img src="../../_static/diagrams/compressor.svg" alt="Compressor diagram" style="max-width:100%">

Driven by a Flownex-style (`.cop`) characteristic map instead of a fixed
`PR`/`eta_s` — pressure ratio and isentropic efficiency are read off the
map's iso-speed curves at the compressor's shaft speed `N` and the corrected
mass flow implied by the current inlet state.

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters:** `map_path`, `gamma`
(optional), `N` (shaft speed `[rev/min]`; leave `None` to solve for it),
`factor_overrides` (optional map calibration), `heat_path` (optional)

$$
A = \frac{N/60}{\sqrt{T_\text{in}}} \quad\text{(corrected speed)}
\qquad
B = \frac{\dot m\sqrt{T_\text{in}}}{P_\text{in}} \quad\text{(corrected flow)}
$$

$$
PR = \text{map.pressure\_ratio}(A, B) \qquad \eta_s = \text{map.efficiency}(A, B)
$$

then the same isentropic relation as [`SimpleCompressor`](simple-compressor.md)
using this `PR` and `eta_s`. **Free parameter:** leave `N=None` to solve for
shaft speed instead of a fixed `PR`/power (needs a matching
[`Setpoint`](../control/setpoint.md)/[`Controller`](../control/controller.md)/
[`Shaft`](../mechanical/shaft.md) residual). `factor_overrides` lets you
calibrate the map's own conversion factors against test data without
editing the map file.

---
Part of [Turbomachinery](index.md).
