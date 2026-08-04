# Heat exchangers & phase change

Single- and two-stream heat exchangers, plus the two-phase (boiling/condensing)
equipment that models heat transfer via saturation enthalpies and quality
instead of `cp` — required once a stream crosses into the two-phase dome,
where `cp → ∞`. Two-phase components require a real-fluid model exposing
saturation properties (`CoolPropFluid`, the `coolprop` extra).

```{toctree}
:maxdepth: 1

simple-heat-exchanger
heat-exchanger
multi-pass-heat-exchanger
condenser
simple-condenser
evaporator
simple-evaporator
drum
tank
```

| Component | Role |
|---|---|
| [`SimpleHeatExchanger`](simple-heat-exchanger.md) | Single-stream, signed duty `Q` (positive heats, negative cools). |
| [`HeatExchanger`](heat-exchanger.md) | Two-stream: fixed effectiveness, or effectiveness derived from `UA` + flow arrangement (auto-selected by which one you give). |
| [`MultiPassHeatExchanger`](multi-pass-heat-exchanger.md) | `HeatExchanger` subclass that's always UA/NTU mode, for when you want that explicit at the call site. |
| [`Condenser`](condenser.md) | Two-stream: working fluid condenses against an explicit coolant stream. |
| [`SimpleCondenser`](simple-condenser.md) | Single-stream: condenses to a spec outlet, coolant not modeled. |
| [`Evaporator`](evaporator.md) | Two-stream: working fluid boils against an explicit heat-source stream. |
| [`SimpleEvaporator`](simple-evaporator.md) | Single-stream: boils to a spec outlet, heat source not modeled. |
| [`Drum`](drum.md) | Steam drum — two-phase reservoir with liquid-level dynamics. |
| [`Tank`](tank.md) | Constant-volume plenum with real mass/energy storage (single-phase). |
