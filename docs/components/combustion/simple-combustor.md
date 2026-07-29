# SimpleCombustor

<img src="../../_static/diagrams/simple_combustor.svg" alt="SimpleCombustor diagram" style="max-width:100%">

A fixed lower-heating-value (LHV) heat-release model: no combustion
chemistry, fuel is treated as pure heat release into the existing working
fluid. Fuel mass is added to the flow.

**Ports:** `in`, `out` (plus `fuel_in` if `use_fuel_port=True`, shown below)
&nbsp;·&nbsp; **Parameters:** `LHV`, `PR` (default 0.97), `efficiency`
(default 1.0), `mdot_fuel` (leave `None` to solve for it), `use_fuel_port`,
`heat_path` (optional)

$$
Q = \dot m_\text{fuel}\cdot LHV \cdot \eta
\qquad
\dot m_\text{out} = \dot m_\text{in} + \dot m_\text{fuel}
$$

$$
P_\text{out} - PR\cdot P_\text{in} = 0
$$

$$
\dot m_\text{out}\,h_\text{out} - \left(\dot m_\text{in}\,h_\text{in} + \dot m_\text{fuel}\,h_\text{fuel,in} + Q - Q_\text{loss}\right) = 0
$$

(`h_fuel,in` only enters if `use_fuel_port=True`; otherwise fuel's own
sensible enthalpy is ignored as negligible next to LHV.)

<img src="../../_static/diagrams/simple_combustor_fuel_port.svg" alt="SimpleCombustor with a fuel port" style="max-width:100%">

**Two ways to get `mdot_fuel`:**

- `use_fuel_port=False` (default): give `mdot_fuel` directly, or leave
  `None` to make it a free parameter driven by a
  [`Setpoint`](../control/setpoint.md)/[`Controller`](../control/controller.md)
  elsewhere (e.g. a target turbine inlet temperature).
- `use_fuel_port=True`: a genuine third port, `fuel_in`, appears — connect a
  real fuel-supply branch to it (`Source → Pipe → combustor`). `mdot_fuel`
  is then read directly from that branch's own solved flow, and its actual
  `(P, h)` is used for a more complete energy balance.

For the same physical role using real combustion chemistry instead of a
fixed LHV, see [`Combustor`](combustor.md).

---
Part of [Combustion](index.md).
