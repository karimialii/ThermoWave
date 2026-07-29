# Sample maps

Ready-to-load example data for the map-driven components
(`Compressor`/`Turbine` via `thermowave.maps.CharacteristicMap`,
`Generator` via `thermowave.maps.TorqueSpeedMap`). Load any of these with
`CharacteristicMap.from_file(...)` / `TorqueSpeedMap.from_file(...)` to try a
component without writing your own map file first; swap in real vendor/rig
data the same way once you have it.

- `compressor_generic.cop` — single constant-geometry compressor map
  (Flownex `.cop` format), one "Angle 0" block per section.
- `turbine_generic.tur` — single constant-geometry turbine map (Flownex
  `.tur` format), same single-angle structure.
- `compressor_vgv_generic.cop` — synthetic variable-geometry compressor map:
  three vane-angle blocks (-10, 0, 10) per section, same file format as
  `compressor_generic.cop` but with `CharacteristicMap` interpolating across
  `angle` as a third dimension (see `CharacteristicMap.pressure_ratio(A, B,
  angle=...)`/`efficiency(A, B, angle=...)`). Synthetic — invented to
  exercise the multi-angle parser/interpolation, not measured rig data.
- `generator_generic.gen` — single speed-vs-torque curve (`TorqueSpeedMap`'s
  plain-text format: `speed torque` per line).

No map-file format exists for `MultiPassHeatExchanger` — its geometry input
is `UA` (plus `arrangement`/`n_passes`, or a `correlation=` callable for a
custom effectiveness relation), not a lookup table, so there's nothing to
sample here for it. See the component's own docstring
(`src/thermowave/components/multi_pass_heat_exchanger.py`).
