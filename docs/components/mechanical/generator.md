# Generator

<img src="../../_static/diagrams/generator.svg" alt="Generator diagram" class="component-diagram">

Two variants, both passive readers with no flow ports and no residuals —
they report a derived electrical output, they don't feed back into the
solve, and neither couples back into the shaft's own power balance (for a
generator that genuinely loads the shaft, a variable-speed machine settling
at whatever speed covers a commanded electrical demand, use
[`ShaftLoad`](shaft-load.md) instead). [`Generator`](#generator-torque-map)
reads torque off a speed-vs-torque characteristic map;
[`SimpleGenerator`](#simplegenerator-fixed-efficiency) applies a fixed
efficiency to shaft power directly.

## `Generator` (torque map)

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

## `SimpleGenerator` (fixed efficiency)

<img src="../../_static/diagrams/simple_generator.svg" alt="SimpleGenerator diagram" class="component-diagram">

Reads power through a `"power"` signal port instead of holding a direct
object reference.

**Ports:** none (flow) — mechanical: `shaft` &nbsp;·&nbsp; signal: `power`
&nbsp;·&nbsp; **Parameters:** `efficiency`

```python
gen = SimpleGenerator(name="gen1", efficiency=0.96)
network.add_component(gen)
network.connect(gen, "power", turbine, "power", kind="signal")
network.connect(gen, "shaft", turbine, "shaft", kind="mechanical")  # optional, for N in report_metrics()
```

$$P_\text{elec} = P_\text{shaft} \cdot \eta$$

---
Part of [Mechanical & electrical](index.md).
