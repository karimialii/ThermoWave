# Control & instrumentation

The components that close free parameters left `None` elsewhere in the
network — measurement, steady targets, time-domain control, and scheduled
profiles. See [The modeling framework](../index.md#free-parameters-giving-a-target-instead-of-a-direct-input)
for how free parameters and their closing residuals fit together.

```{toctree}
:maxdepth: 1

sensor
controller
pid-controller
setpoint
schedule
```

| Component | Role |
|---|---|
| [`Sensor`](sensor.md) | A passive measurement tap on a network node. |
| [`Controller`](controller.md) | Drives a free parameter until a `Sensor` reads a target (ideal, infinite-gain). |
| [`PIDController`](pid-controller.md) | Time-domain PID counterpart to `Controller`, for `solve_transient()`. |
| [`Setpoint`](setpoint.md) | Drives a free parameter until the *same* component's own metric hits a target. |
| [`Schedule`](schedule.md) | Drives a plain attribute along a piecewise time profile during `solve_transient()`. |
