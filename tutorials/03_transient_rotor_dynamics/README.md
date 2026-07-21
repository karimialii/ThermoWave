# Tutorial 3: Transient rotor dynamics

Gives the shaft real inertia and watches it spool up over time, instead of
assuming it's already at equilibrium.

| Step | File | Teaches |
|---|---|---|
| 1 | `01_dynamic_shaft_equilibrium.py` | `Shaft(dynamic=True, inertia=...)` makes speed a differential state; a steady solve finds the torque-balance equilibrium with no Setpoint needed. |
| 2 | `02_spool_up_transient.py` | `Network.solve_transient()` integrates that state forward from an off-equilibrium start, producing a genuine spool-up curve (plotted). |

Run any step directly, e.g.:

```
.venv/bin/python tutorials/03_transient_rotor_dynamics/02_spool_up_transient.py
```

Prerequisite: [Tutorial 2](../02_gas_turbine_cycle/) (shaft coupling, in its
steady-state/Setpoint-driven form).

**Next:** [Tutorial 4: Closed-loop control](../04_closed_loop_control/) — a
different kind of time-varying behavior, where PID setpoints change rather
than the shaft's own physics.
