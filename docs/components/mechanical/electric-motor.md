# ElectricMotor

<img src="../../_static/diagrams/electric_motor.svg" alt="ElectricMotor diagram" style="max-width:100%">

Electrically-driven mechanical power source — the inverse of
[`SimpleGenerator`](simple-generator.md). A passive reader with no flow
ports: it reads a mechanical component's own required
`report_metrics()["power [W]"]` (e.g. an electrically-driven compressor with
no turbine on its shaft) and reports the electrical power that must be
drawn to supply it.

**Ports:** none &nbsp;·&nbsp; **Parameters:** `component` (the mechanical
load), `efficiency`

$$
P_\text{elec} = \frac{P_\text{shaft, required}}{\eta}
$$

**Residuals:** none — purely a derived reading. The mechanical component's
own free speed, if it has one, still needs its own
[`Setpoint`](../control/setpoint.md)/[`Controller`](../control/controller.md).

---
Part of [Mechanical & electrical](index.md).
