# Changelog

Notable changes to ThermoWave, in reverse chronological order. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions
correspond to the `[project] version` in `pyproject.toml`.

## 0.8.0 — 2026-08-04

### Added

- **Mechanical and signal connection kinds.** `network.connect()` now
  supports `kind="mechanical"` (shares a shaft speed `N`) and
  `kind="signal"` (shares a scalar one component computes and another
  reads, e.g. power) alongside `kind="flow"`, resolved through the same
  union-find as flow ports (`Network._all_mechanical_nodes()`,
  `NetworkState.N()`/`.signal()`). `Turbine`/`Compressor` expose a
  `"shaft"` mechanical port and a `"power"` signal port; `Generator`/
  `ElectricMotor`/`SimpleGenerator` read speed/power through
  `connect()`-wired ports instead of a direct object reference.
  `Setpoint`/`Controller`/`PIDController` accept a mechanical port name
  (e.g. `free_param="shaft"`) anywhere they used to accept only an
  ordinary free parameter. See
  [Mechanical & electrical](components/mechanical/index.md) and
  [the gas-turbine tutorial](tutorials/building-a-gas-turbine-model.md).
- **`HeatExchanger`**, a new two-stream heat exchanger that auto-selects
  fixed-effectiveness mode (`effectiveness=...`) or UA/NTU/arrangement mode
  (`UA=...`) by which one you give. See
  [`HeatExchanger`](components/heat-exchangers/heat-exchanger.md).
- **`Combustor(use_fuel_port=True)`** adds a genuine `fuel_in` port
  (mirroring `SimpleCombustor`'s existing flag) — `mdot_fuel` is read from
  that port's own solved mass flow, and if the connected fuel stream's
  fluid is itself Cantera-flavored, its actual composition and its own
  `(P, h)` are used in the equilibrium calculation instead of the `fuel=`
  string evaluated at the air inlet's state. See
  [`Combustor`](components/combustion/combustor.md).

### Changed

- **`SimpleHeatExchanger` is now single-stream.** It has 2 ports (`in`,
  `out`) driven by a signed duty `Q` [W] (positive heats, negative cools)
  instead of the old 4-port (`hot_in`/`hot_out`/`cold_in`/`cold_out`)
  fixed-effectiveness dual-fluid model. **Migration:** replace
  `SimpleHeatExchanger(effectiveness=..., PR_hot=..., PR_cold=...)` with
  `HeatExchanger(effectiveness=..., PR_hot=..., PR_cold=...)`. See
  [`SimpleHeatExchanger`](components/heat-exchangers/simple-heat-exchanger.md).
- **`MultiPassHeatExchanger` moved** from
  `thermowave.components.multi_pass_heat_exchanger` to
  `thermowave.components.heat_exchanger`, where it's now a thin subclass of
  `HeatExchanger` (constructor and behavior unchanged — no math is
  duplicated between the two classes). The flat `thermowave.components`
  import (`from thermowave.components import MultiPassHeatExchanger`)
  is unaffected.

### Breaking changes

- **`Shaft`'s constructor**: `components=` is renamed `members=`, and it no
  longer reads members' free parameters at construction time. Every
  mechanical/power coupling now needs an explicit
  `network.connect(shaft, "m{i}", member, "shaft", kind="mechanical")` /
  `network.connect(shaft, "p{i}", member, "power", kind="signal")` call
  after `add_component()` — a member is no longer speed-tied just by
  appearing in the list. See
  [`Shaft`](components/mechanical/shaft.md) for the updated wiring pattern.
- **`Generator`, `ElectricMotor`, `SimpleGenerator`**: the `component=`
  constructor argument is gone. They now read speed/power through
  `connect()`-wired `"shaft"` (mechanical) and `"power"` (signal) ports
  instead of holding a direct object reference — wire them the same way as
  any other mechanical/signal port. See
  [Mechanical & electrical](components/mechanical/index.md).
- **`Turbine`/`Compressor`**: `N` is no longer reported via
  `free_parameters()`/read via `state.param(f"{name}.N")` — it's now a
  `"shaft"` mechanical port (fixed when `N` is given, a free Newton
  unknown when `N=None`), read via `state.N(...)`. Any code calling
  `component.free_parameters()` expecting an `"N"` key, or reading
  `result.params[f"{name}.N"]` directly, should instead check
  `component.free_mechanical_ports()` / `component.mechanical_ports()`
  and read `result.node_N[f"{name}.shaft"]`.
- **`SimpleHeatExchanger`**: see *Changed*, above — this is also a
  breaking API/behavior change, not just a rename.

---

Earlier releases aren't yet broken out into dated entries here — see the
[README's Roadmap section](https://github.com/karimialii/ThermoWave#roadmap)
and `git log` for the fuller history up to this point.
