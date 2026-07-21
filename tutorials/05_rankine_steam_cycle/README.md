# Tutorial 5: Rankine steam cycle

Builds a complete steam power cycle on real water/steam properties
(`CoolPropFluid`), one component at a time.

| Step | File | Teaches |
|---|---|---|
| 1 | `01_pump_and_boiler.py` | `Pump` (entropy-based pressurization) and `SimpleEvaporator` (boil + superheat to a target outlet condition). |
| 2 | `02_add_steam_turbine.py` | `SteamTurbine` — the wet-steam-correct expansion, reporting exhaust quality. |
| 3 | `03_full_cycle_with_condenser.py` | `SimpleCondenser` closes the loop; computing cycle efficiency from the four corner enthalpies. |

Run any step directly, e.g.:

```
.venv/bin/python tutorials/05_rankine_steam_cycle/03_full_cycle_with_condenser.py
```

No prerequisite tutorial — this series is self-contained (it doesn't build
on the gas-turbine tutorials). Requires the `coolprop` extra:
`pip install thermowave[coolprop]`.

**Next:** [Tutorial 6: Boiler drum transient](../06_boiler_drum_transient/) —
the steam drum that would sit upstream of this cycle's boiler in a real
plant, and the level dynamics that come with it.
