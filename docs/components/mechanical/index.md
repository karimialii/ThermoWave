# Mechanical & electrical

The components that couple turbomachinery together mechanically, and that
convert mechanical shaft power to and from electrical power. None of these
have flow ports — they read other components' `report_metrics()` and/or
contribute speed-tie residuals, rather than moving fluid. See
[Turbomachinery](../turbomachinery/index.md) for the compressors/turbines/
pumps these attach to.

```{toctree}
:maxdepth: 1

shaft
shaft-load
electric-motor
generator
simple-generator
```

| Component | Role |
|---|---|
| [`Shaft`](shaft.md) | Ties two or more machines' shaft speeds together (steady or dynamic/rotor-inertia). |
| [`ShaftLoad`](shaft-load.md) | A genuine torque on a `Shaft` — the electrical demand a generator places on the rotor. |
| [`ElectricMotor`](electric-motor.md) | Reports the electrical power drawn to supply an electrically-driven mechanical load. |
| [`Generator`](generator.md) | Passive reader: shaft speed → torque (map) → electrical power. |
| [`SimpleGenerator`](simple-generator.md) | Passive reader: shaft power × fixed efficiency → electrical power. |
