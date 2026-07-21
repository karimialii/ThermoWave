# Tutorial 2: Gas turbine cycle

Builds a recuperated, single-shaft gas turbine one concept at a time, using
the same T100 microturbine compressor/turbine maps throughout.

| Step | File | Teaches |
|---|---|---|
| 1 | `01_basic_brayton_cycle.py` | Compressor -> combustor -> Turbine with independent, fixed speeds — the basic thermodynamic cycle. |
| 2 | `02_shaft_coupling.py` | `Shaft` ties compressor and turbine to one physical speed; `Setpoint` closes the one remaining degree of freedom. |
| 3 | `03_recuperated_cycle.py` | A recuperator recovers exhaust heat, and fuel arrives through a real physical supply line instead of an abstract number. |

Run any step directly, e.g.:

```
.venv/bin/python tutorials/02_gas_turbine_cycle/03_recuperated_cycle.py
```

Prerequisite: [Tutorial 1](../01_building_a_flow_network/) (network basics).

**Next:**
- [Tutorial 3: Transient rotor dynamics](../03_transient_rotor_dynamics/) — the
  shaft's own inertia, spooling up over time instead of sitting at
  equilibrium.
- [Tutorial 4: Closed-loop control](../04_closed_loop_control/) — replacing
  this cycle's fixed setpoints with real PID control loops.
