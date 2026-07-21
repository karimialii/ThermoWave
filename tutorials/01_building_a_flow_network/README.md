# Tutorial 1: Building a flow network

The core mechanics every ThermoWave network uses, regardless of how complex
it eventually gets: pick a fluid, add components, connect their ports,
solve.

| Step | File | Teaches |
|---|---|---|
| 1 | `01_source_and_sink.py` | The minimum network: `Network(fluid=...)`, `add_component()`, `connect()`, `solve()`. |
| 2 | `02_full_chain.py` | Chaining components (compressor, pipe, valve) — the pattern doesn't change as the network grows. |

Run either script directly:

```
.venv/bin/python tutorials/01_building_a_flow_network/01_source_and_sink.py
.venv/bin/python tutorials/01_building_a_flow_network/02_full_chain.py
```

**Next:** [Tutorial 2: Gas turbine cycle](../02_gas_turbine_cycle/) — a
turbomachinery cycle with a compressor and turbine sharing one shaft.
