# Tutorials

Step-by-step guides for building complex systems with ThermoWave, each
broken into a short sequence of runnable scripts that add one new concept
at a time. For a one-file-per-component reference instead, see
[`examples/`](../examples/).

## Suggested order

1. **[01_building_a_flow_network](01_building_a_flow_network/)** — the core
   mechanics every network uses: `Network`, `add_component()`, `connect()`,
   `solve()`. Start here if you're new to ThermoWave.
2. **[02_gas_turbine_cycle](02_gas_turbine_cycle/)** — a recuperated,
   single-shaft gas turbine: map-based turbomachinery, `Shaft` coupling, a
   real fuel supply line.
3. **[03_transient_rotor_dynamics](03_transient_rotor_dynamics/)** — giving
   that shaft real inertia and watching it spool up over time.
4. **[04_closed_loop_control](04_closed_loop_control/)** — replacing fixed
   setpoints with PID control loops, then a time-varying setpoint profile.
5. **[05_rankine_steam_cycle](05_rankine_steam_cycle/)** — a complete steam
   power cycle on real water/steam properties (two-phase components).
   Self-contained; doesn't require 1-4.
6. **[06_boiler_drum_transient](06_boiler_drum_transient/)** — the steam
   drum that would sit upstream of tutorial 5's boiler, and its level
   dynamics under a load step.

Tutorials 2-4 build on each other (one continuous gas-turbine cycle,
progressively extended); 5-6 are a separate, self-contained steam/two-phase
track. Either track can be started directly from tutorial 1.

Each tutorial folder has its own `README.md` with a per-step breakdown.
Every script is directly runnable and prints what to look at; most also
tell you what the next step adds.

Some tutorials require an optional extra:
- Tutorials 5-6 need `pip install thermowave[coolprop]` (real water/steam
  properties).
- Plotting (`history.plot(...)`) needs `pip install thermowave[plot]`,
  gracefully skipped where it isn't installed.
