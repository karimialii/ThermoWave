# SimpleCompressor

<img src="../../_static/diagrams/simple_compressor.svg" alt="SimpleCompressor diagram" style="max-width:100%">

Analytic compressor: fixed pressure ratio and isentropic efficiency, no
performance map.

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters:** `PR`, `eta_s`, `gamma`
(optional, else derived from the fluid), `heat_path` (optional)

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
fully adiabatic (`Q_loss = 0`).

For the same physics driven by a real characteristic map instead of a fixed
`PR`/`eta_s`, see [`Compressor`](compressor.md).

---
Part of [Turbomachinery](index.md).
