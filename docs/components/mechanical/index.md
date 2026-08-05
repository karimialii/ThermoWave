# Mechanical & electrical

The components that couple turbomachinery together mechanically, and that
convert mechanical shaft power to and from electrical power. None of these
have flow ports — they couple through **mechanical** (`kind="mechanical"`,
a shared shaft speed `N`) and **signal** (`kind="signal"`, a scalar one
component computes and another reads, e.g. power) ports instead, wired with
`network.connect(a, "shaft", b, "shaft", kind="mechanical")` /
`network.connect(a, "power", b, "power", kind="signal")` the same way flow
ports are, rather than holding direct object references. See
[Turbomachinery](../turbomachinery/index.md) for the compressors/turbines/
pumps these attach to (`Turbine`/`Compressor` expose their own `"shaft"`
mechanical port and `"power"` signal port).

When the mechanical connection goes through a [`Shaft`](shaft.md)
specifically, the matching signal connection is made automatically the
moment the mechanical one is — see `Shaft.on_connected()`. The two-call
form above is still what you need for a direct machine-to-machine
mechanical link with no `Shaft` in between (e.g. `Generator`/
`ElectricMotor` reading straight off a turbine's own `"shaft"` port).

```{toctree}
:maxdepth: 1

shaft
shaft-load
electric-motor
generator
```

| Component | Role |
|---|---|
| [`Shaft`](shaft.md) | Ties two or more machines' shaft speeds together (steady or dynamic/rotor-inertia). |
| [`ShaftLoad`](shaft-load.md) | A genuine torque on a `Shaft` — the electrical demand a generator places on the rotor. |
| [`ElectricMotor`](electric-motor.md) | Reports the electrical power drawn to supply an electrically-driven mechanical load. |
| [`Generator`](generator.md#generator-torque-map) | Passive reader: shaft speed → torque (map) → electrical power. |
| [`SimpleGenerator`](generator.md#simplegenerator-fixed-efficiency) | Passive reader: shaft power × fixed efficiency → electrical power. |
