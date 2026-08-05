# Turbomachinery

The fluid-side machines — compressors, turbines, and pumps — that add or
extract work from the flow. For the shafts, loads, motors, and generators
that couple these machines mechanically/electrically, see
[Mechanical & electrical](../mechanical/index.md).

```{toctree}
:maxdepth: 1

compressor
turbine
pump
```

| Component | Role |
|---|---|
| [`Compressor`](compressor.md#compressor-map-driven) | Map-driven (Flownex `.cop`), `PR`/efficiency read off iso-speed curves. |
| [`SimpleCompressor`](compressor.md#simplecompressor-fixed-ratio) | Fixed pressure ratio + isentropic efficiency, no map. |
| [`Turbine`](turbine.md#turbine-map-driven) | Map-driven (Flownex `.tur`), same corrected-parameter approach as `Compressor`. |
| [`SimpleTurbine`](turbine.md#simpleturbine-fixed-ratio) | Fixed expansion ratio + isentropic efficiency, no map. |
| [`SteamTurbine`](turbine.md#steamturbine-wet-steam) | Entropy-based expansion, correct through the two-phase dome. |
| [`Pump`](pump.md) | Entropy-based liquid pump, closes a Rankine loop's low-pressure side. |
