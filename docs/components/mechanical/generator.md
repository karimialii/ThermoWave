# Generator

<img src="../../_static/diagrams/generator.svg" alt="Generator diagram" style="max-width:100%">

A passive reader of a shaft component's own `report_metrics()["N [rev/min]"]`,
with no flow ports and no residuals — it reports a derived electrical
output, it doesn't feed back into the solve. Driven by a speed-vs-torque
characteristic map (`TorqueSpeedMap`, `.gen`-style file) instead of a fixed
efficiency:

**Ports:** none &nbsp;·&nbsp; **Parameters:** `component` (the shaft to
read speed from), `map_path`, `efficiency` (default 1.0)

$$
\omega = \frac{2\pi N}{60}
\qquad
P_\text{mech} = \text{map.torque}(N)\cdot\omega
\qquad
P_\text{elec} = P_\text{mech}\cdot\eta
$$

Does not couple back into the shaft's own power balance — for a generator
that genuinely loads the shaft (a variable-speed machine settling at
whatever speed covers a commanded electrical demand), use
[`ShaftLoad`](shaft-load.md) instead. For a fixed-efficiency reader instead
of a torque map, see [`SimpleGenerator`](simple-generator.md).

---
Part of [Mechanical & electrical](index.md).
