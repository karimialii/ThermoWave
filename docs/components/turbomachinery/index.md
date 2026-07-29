# Turbomachinery

The fluid-side machines — compressors, turbines, and pumps — that add or
extract work from the flow. For the shafts, loads, motors, and generators
that couple these machines mechanically/electrically, see
[Mechanical & electrical](../mechanical/index.md).

```{toctree}
:maxdepth: 1

simple-compressor
compressor
simple-turbine
turbine
steam-turbine
pump
```

| Component | Role |
|---|---|
| [`SimpleCompressor`](simple-compressor.md) | Fixed pressure ratio + isentropic efficiency, no map. |
| [`Compressor`](compressor.md) | Map-driven (Flownex `.cop`), `PR`/efficiency read off iso-speed curves. |
| [`SimpleTurbine`](simple-turbine.md) | Fixed expansion ratio + isentropic efficiency, no map. |
| [`Turbine`](turbine.md) | Map-driven (Flownex `.tur`), same corrected-parameter approach as `Compressor`. |
| [`SteamTurbine`](steam-turbine.md) | Entropy-based expansion, correct through the two-phase dome. |
| [`Pump`](pump.md) | Entropy-based liquid pump, closes a Rankine loop's low-pressure side. |
