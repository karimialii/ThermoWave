# SimpleTurbine

<img src="../../_static/diagrams/simple_turbine.svg" alt="SimpleTurbine diagram" style="max-width:100%">

Analytic turbine: fixed expansion ratio (`PR = P_in/P_out > 1`) and
isentropic efficiency.

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters:** `PR`, `eta_s`, `gamma`
(optional), `heat_path` (optional)

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

For the same physics driven by a real characteristic map, see
[`Turbine`](turbine.md); for expansion through the two-phase dome (wet
steam), see [`SteamTurbine`](steam-turbine.md).

---
Part of [Turbomachinery](index.md).
