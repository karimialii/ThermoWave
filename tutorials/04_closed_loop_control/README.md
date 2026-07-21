# Tutorial 4: Closed-loop control

Replaces the exact, instantaneous `Setpoint` from tutorial 2 with real
time-domain PID control loops, building up to several loops plus a
time-varying setpoint driven by a `Schedule`.

| Step | File | Teaches |
|---|---|---|
| 1 | `01_single_pid_loop.py` | One `PIDController` loop (fuel flow -> turbine outlet temperature) settling toward a target over simulated time. |
| 2 | `02_two_pid_loops.py` | A second, independent loop (shaft speed -> electrical power) running alongside the first — and why a loop's gain sign follows the system's actual sensitivity, not intuition. |
| 3 | `03_power_profile_schedule.py` | `Schedule` drives one loop's setpoint through a hold/ramp/hold profile, still inside a single `solve_transient()` call. |

Run any step directly, e.g.:

```
.venv/bin/python tutorials/04_closed_loop_control/03_power_profile_schedule.py
```

Prerequisite: [Tutorial 2](../02_gas_turbine_cycle/) (the cycle these loops
control) and [Tutorial 3](../03_transient_rotor_dynamics/) (`solve_transient()`
basics).

This tutorial covers the same ground as the T100 microturbine power-profile
example that originally lived in `examples/`, using `SimpleCombustor`
(fixed lower heating value) instead of the Cantera chemical-equilibrium
`Combustor` — no optional `cantera` dependency required, and no need for the
staged-continuation bootstrapping a harder multi-target cold start needs.
