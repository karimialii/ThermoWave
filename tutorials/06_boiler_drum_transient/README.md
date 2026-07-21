# Tutorial 6: Boiler drum transient

Introduces the `Drum`'s differential-state mechanism in isolation, then
integrates it through a real disturbance.

| Step | File | Teaches |
|---|---|---|
| 1 | `01_drum_state_derivative.py` | Calling `state_derivative()` directly at hand-picked balanced/draining states, no network or time-stepping needed. |
| 2 | `02_level_response_to_demand_step.py` | Integrating that state forward through a steam-demand step, producing a real level/pressure trajectory (plotted). |

Run any step directly, e.g.:

```
.venv/bin/python tutorials/06_boiler_drum_transient/02_level_response_to_demand_step.py
```

Prerequisite: none required, though [Tutorial 3](../03_transient_rotor_dynamics/)
covers the same differential-state mechanism (`differential_parameters()`/
`state_derivative()`) on a `Shaft` first, which may make this one land
faster. Requires the `coolprop` extra: `pip install thermowave[coolprop]`.

A `Drum`'s level has no steady-state value — it's a pure integrator, the
same reason `Network.solve()` is singular in a drum's own level (see
`Drum`'s own docstring). That's why this tutorial works directly with
`state_derivative()` and manual time integration rather than
`Network.solve_transient()` on a full network, which would need a level
controller (see [Tutorial 4](../04_closed_loop_control/)) to be well-posed.
