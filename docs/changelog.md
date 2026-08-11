# Changelog

Notable changes to ThermoWave, in reverse chronological order. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions
correspond to the `[project] version` in `pyproject.toml`.

## 0.9.2 — 2026-08-11

### Docs

- Every benchmark's runnable script now lives in its own standalone public
  repo (`karimialii/<name>-benchmark`) instead of a gitignored copy at the
  package repo root — each is independently cloneable, with its own
  README, license, and plotting code. The docs pages here keep the full
  written analysis and pre-rendered plots, linking out to the
  corresponding repo for "run it yourself."
- New component diagrams for `Recycle` and `SimpleHeater` (previously
  text-only pages); refreshed `FeedwaterHeater`/`SimpleCondenser`/
  `SimpleEvaporator` artwork; the SEGS and ORC benchmark pages now embed
  their own plant/cycle layout diagrams.
- No source code changes in this release — version bump is so Read the
  Docs' `stable` alias (which tracks the latest tag) picks up the above.

## 0.9.1 — 2026-08-11

### Fixed

- **Exergy analysis silently skipped `FeedwaterHeater`.** `FeedwaterHeater`
  is a genuine two-stream exchanger (it shares `HeatExchanger`'s default
  `hot_in`/`hot_out`/`cold_in`/`cold_out` port names) but wasn't registered
  in `thermowave.core.exergy`'s two-stream dispatch, so it fell through
  every branch with no `E_F`/`E_P`/`E_D` row and no contribution to the
  network totals — the same failure mode `Condenser` needed fixing for
  previously. Now costed identically to any other two-stream exchanger.

### Docs

- SEGS benchmark doc/plots no longer reference TESPy — every design value
  (`ttd_u`/`ttd_l`) is now framed as a standard closed-FWH design target,
  not a TESPy-sourced one. Adds the plant layout diagram and two new
  plots (solved FWH `ttd` vs. target, solved vs. Table 3 extraction
  split).

## 0.9.0 — 2026-08-10

### Added

- **`HumidAirFluid`**, a psychrometric fluid model backed by CoolProp's
  `HAPropsSI` — the first fluid model in ThermoWave that tracks humidity
  ratio, relative humidity, wet-bulb temperature, and the latent/sensible
  energy split. Humidity ratio `W` [kg water/kg dry air] is fixed at
  construction (the same role `CanteraFluid`'s composition plays); `mdot`
  for a `HumidAirFluid` stream means **dry-air** mass flow, the ASHRAE
  convention, and the basis `HAPropsSI`'s own `H`/`Vha`/`cp_ha` outputs
  already use. `from_relative_humidity(T, RH, P)` derives `W` from a more
  common "50% RH at 20°C" specification. Drops into any existing
  `Source(fluid=...)`/`Pipe`/`HeatExchanger`/`Junction` network with no
  other changes required. See `thermowave.fluids.humid_air.HumidAirFluid`
  and `thermowave.fluids.psychrometrics.supports_humid_air()`.
- **`CoolingTower`**, evaporative cooling with real water-vapor mass
  transfer — a water stream loses heat *and* mass, an air stream gains
  heat *and* humidity, closed by a saturation constraint on the outlet air
  (`target_RH_out`, default 1.0). See
  [`CoolingTower`](components/heat-exchangers/cooling-tower.md).
- **`FeedwaterHeater`** and **`Deaerator`**, first-class regenerative
  feedwater heater components — replacing the pattern of hand-building each
  one from a `SimpleCondenser` + `SimpleEvaporator` + `Junction` triple with
  two independently-fixed (and never exactly matching) duties.
  `FeedwaterHeater` (closed FWH) is structurally `Condenser` renamed and
  generalized — one shared duty computed once, not two. `Deaerator` (open
  FWH) is `Junction`'s own mixing plus the one constraint `Junction` lacks:
  pinning the outlet to saturated liquid. See
  [Feedwater heaters](components/heat-exchangers/feedwater-heaters.md).

### Changed

- **mypy noise on fluid-capability duck typing eliminated.** The
  `two_phase.py`/`psychrometrics.py` runtime capability checks
  (`supports_two_phase()`, `supports_entropy()`, `supports_humid_air()`, and
  the new `supports_cantera_composition()`) are now typed as `TypeIs`
  Protocols, removing ~40 recurring false-positive "`BaseFluid` has no
  attribute" errors. Purely an internal type-checking improvement — no
  runtime behavior change. Adds `typing_extensions` as a direct dependency
  (needed for `TypeIs` on the package's Python 3.10 floor).

### Fixed

- **`SteamTurbine.report_metrics()`** could raise a bare `AttributeError`
  for a fluid that supports entropy but not two-phase queries (e.g. a
  hypothetical future ideal-gas-only fluid) — now raises the same clear
  `require_two_phase()` error every other two-phase call site already has.

## 0.8.2 — 2026-08-05

### Added

- **`Shaft` wires power automatically.** Connecting a member's mechanical
  `"shaft"` port (`network.connect(shaft, "m{i}", member, "shaft",
  kind="mechanical")`) now also wires that same member's `"power"` signal
  port to the matching `"p{i}"` — the new `BaseComponent.on_connected()`
  hook, fired by `Network.connect()` on both endpoints of every
  connection, with `Shaft` as its only override. A torque-only member
  (`ShaftLoad`, no `"shaft"` mechanical port) still needs its own explicit
  signal `connect()` call, since there's no mechanical connection to
  piggyback on. See [`Shaft`](components/mechanical/shaft.md).
- **`Shaft.autowire(network)`**, called once after `add_component(shaft)`,
  makes every mechanical/signal `connect()` call a shaft's members need in
  one go — for any member whose ports are named the usual `"shaft"`/
  `"power"` way (every built-in rotating/electrical component already is).
- **`BaseComponent.shaft_sign()`** — a fixed per-class default for
  `Shaft`'s `signs`: `+1.0` (delivers power) on `Turbine`/`ElectricMotor`,
  `-1.0` (draws power) on `Compressor`/`Generator`/`SimpleGenerator`/
  `ShaftLoad`. `Shaft(signs=None)` (the new default) infers `signs` from
  each member's `shaft_sign()` instead of requiring it spelled out by
  hand; pass `signs=[...]` explicitly to override for a genuine edge case
  (a motor running in regenerative mode, a back-driven compressor).
  `Shaft` raises a clear error naming any member with no default
  (`shaft_sign()` returns `None`) rather than silently guessing.

### Changed

- **Read the Docs redesign.** Component diagrams are now centered and
  rendered at a consistent visual size (`docs/_static/custom.css`,
  `object-fit: contain` on a fixed box — a plain `max-width` cap left
  small-intrinsic-size SVGs tiny and everything left-aligned). The
  "Simple"/full component page pairs are merged into one page per family
  (`Compressor`+`SimpleCompressor`, `Turbine`+`SimpleTurbine`+
  `SteamTurbine`, `HeatExchanger`+`SimpleHeatExchanger`+
  `MultiPassHeatExchanger`, `Condenser`+`SimpleCondenser`,
  `Evaporator`+`SimpleEvaporator`, `Combustor`+`SimpleCombustor`,
  `Generator`+`SimpleGenerator`), each as clearly labeled sections on one
  page instead of a separate page per variant — every navigation link and
  cross-reference across the docs was updated to the new page/anchor
  layout.

## 0.8.1 — 2026-08-04

### Added

- **Contributor infrastructure.** `CONTRIBUTING.md` (dev setup, the
  `pytest`/`ruff`/`mypy` commands CI runs, PR conventions),
  `CODE_OF_CONDUCT.md` (Contributor Covenant), and `SECURITY.md` (private
  vulnerability reporting), each mirrored on Read the Docs via
  [Contributing](contributing.md) and [Code of Conduct](code_of_conduct.md).
- CI/Docs/License/Python badges on the README.

### Fixed

- README "Future development" roadmap no longer lists CI and
  `[project.urls]` as outstanding — both landed with 0.8.0's CI workflow
  and haven't been reflected in the doc until now.

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
  [`SimpleHeatExchanger`](components/heat-exchangers/heat-exchanger.md#simpleheatexchanger-single-stream).
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
