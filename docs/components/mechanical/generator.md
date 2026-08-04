# Generator

<img src="../../_static/diagrams/generator.svg" alt="Generator diagram" style="max-width:100%">

A passive reader of shaft speed, with no flow ports and no residuals — it
reports a derived electrical output, it doesn't feed back into the solve.
Reads speed through a real `"shaft"` mechanical port instead of holding a
direct object reference. Driven by a speed-vs-torque characteristic map
(`TorqueSpeedMap`, `.gen`-style file) instead of a fixed efficiency:

**Ports:** none (flow) — mechanical: `shaft` &nbsp;·&nbsp; signal: `power`
&nbsp;·&nbsp; **Parameters:** `map_path`, `efficiency` (default 1.0)

```python
gen = Generator(name="gen1", map_path="alternator.gen")
network.add_component(gen)
network.connect(gen, "shaft", turbine, "shaft", kind="mechanical")
```

$$
\omega = \frac{2\pi N}{60}
\qquad
P_\text{mech} = \text{map.torque}(N)\cdot\omega
\qquad
P_\text{elec} = P_\text{mech}\cdot\eta
$$

Also exposes a `"power"` signal port (`provided_signal_values()` publishes
the same `P_elec` as `report_metrics()["power [W]"]`) for anything
downstream that wants to read it via `kind="signal"` instead of
`report_metrics()`.

Does not couple back into the shaft's own power balance — for a generator
that genuinely loads the shaft (a variable-speed machine settling at
whatever speed covers a commanded electrical demand), use
[`ShaftLoad`](shaft-load.md) instead. For a fixed-efficiency reader instead
of a torque map, see [`SimpleGenerator`](simple-generator.md).

---
Part of [Mechanical & electrical](index.md).
