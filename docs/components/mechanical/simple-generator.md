# SimpleGenerator

<img src="../../_static/diagrams/simple_generator.svg" alt="SimpleGenerator diagram" style="max-width:100%">

A passive reader of a shaft component's own `report_metrics()["power [W]"]`,
with no flow ports and no residuals — it reports a derived electrical
output, it doesn't feed back into the solve.

**Ports:** none &nbsp;·&nbsp; **Parameters:** `component` (the shaft to
read power from), `efficiency`

$$P_\text{elec} = P_\text{shaft} \cdot \eta$$

For torque read off a speed-vs-torque characteristic map instead of a fixed
efficiency, see [`Generator`](generator.md). Neither couples back into the
shaft's own power balance — for a generator that genuinely loads the shaft,
use [`ShaftLoad`](shaft-load.md) instead.

---
Part of [Mechanical & electrical](index.md).
