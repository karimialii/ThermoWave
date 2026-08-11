# Benchmarks

Each benchmark below builds a full ThermoWave network for a real,
published system and checks the solved result against that publication's
own numbers — not a synthetic test case. Every page includes plots of the
solved result against the published reference data, plus the reasoning
and derivations behind the comparison.

The runnable `.py` script behind each page lives in its own standalone
repo (linked from that page, alongside a `git clone` command) rather than
in this one — each is a small, independently cloneable project with its
own README and plotting code.

```{toctree}
:maxdepth: 1

sco2_recompression/sco2_recompression_benchmark
segs_exergy/segs_exergy_benchmark
orc11kwe/orc11kwe_benchmark
rotor37/rotor37_benchmark
hanoi_network/hanoi_trunk_benchmark
```

| Benchmark | System | Reference |
|---|---|---|
| [Recompression sCO2 Brayton cycle](sco2_recompression/sco2_recompression_benchmark.md) | Recuperated supercritical-CO2 power cycle, 5-component exergy analysis | Penkuhn & Tsatsaronis (2018) |
| [30 MWe SEGS solar Rankine cycle](segs_exergy/segs_exergy_benchmark.md) | 7-stage steam turbine, 5 regenerative feedwater heaters, full drain cascade | Lippke, SAND95-1293 (Sandia National Laboratories, 1995) |
| [11 kWe ORC waste-heat rig](orc11kwe/orc11kwe_benchmark.md) | Organic Rankine cycle test rig, R245fa | Published, instrumented test data |
| [NASA Rotor 37](rotor37/rotor37_benchmark.md) | Transonic axial compressor rotor design point | NASA CFD validation case |
| [Hanoi water-distribution network](hanoi_network/hanoi_trunk_benchmark.md) | Trunk-main pipe network | Fujiwara & Khang (1990) |
