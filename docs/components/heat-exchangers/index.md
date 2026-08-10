# Heat exchangers & phase change

Single- and two-stream heat exchangers, plus the two-phase (boiling/condensing)
equipment that models heat transfer via saturation enthalpies and quality
instead of `cp` — required once a stream crosses into the two-phase dome,
where `cp → ∞`. Two-phase components require a real-fluid model exposing
saturation properties (`CoolPropFluid`, the `coolprop` extra).

```{toctree}
:maxdepth: 1

heat-exchanger
condenser
evaporator
heater
drum
tank
feedwater-heaters
cooling-tower
```

| Component | Role |
|---|---|
| [`HeatExchanger`](heat-exchanger.md#heatexchanger-two-stream) | Two-stream: fixed effectiveness, or effectiveness derived from `UA` + flow arrangement (auto-selected by which one you give; leave both unset to also solve for `UA`). |
| [`SimpleHeatExchanger`](heat-exchanger.md#simpleheatexchanger-single-stream) | Single-stream, signed duty `Q` (positive heats, negative cools). |
| [`MultiPassHeatExchanger`](heat-exchanger.md#multipassheatexchanger) | `HeatExchanger` subclass that's always UA/NTU mode, for when you want that explicit at the call site. |
| [`Condenser`](condenser.md#condenser-two-stream) | Two-stream: working fluid condenses against an explicit coolant stream. |
| [`SimpleCondenser`](condenser.md#simplecondenser-single-stream) | Single-stream: condenses to a spec outlet, coolant not modeled. |
| [`Evaporator`](evaporator.md#evaporator-two-stream) | Two-stream: working fluid boils against an explicit heat-source stream. |
| [`SimpleEvaporator`](evaporator.md#simpleevaporator-single-stream) | Single-stream: boils to a spec outlet, heat source not modeled. |
| [`SimpleHeater`](heater.md) | Single-stream: heats to a target *temperature* (not saturation-relative), safe for supercritical fluids. |
| [`Drum`](drum.md) | Steam drum — two-phase reservoir with liquid-level dynamics. |
| [`Tank`](tank.md) | Constant-volume plenum with real mass/energy storage (single-phase). |
| [`FeedwaterHeater`](feedwater-heaters.md#feedwaterheater-closed) | Closed regenerative FWH — condensing extraction steam heats feedwater with one shared duty. |
| [`Deaerator`](feedwater-heaters.md#deaerator-open) | Open regenerative FWH — N-inlet mixing pinned to a saturated-liquid outlet. |
| [`CoolingTower`](cooling-tower.md) | Evaporative cooling — real water-vapor mass transfer to a humid-air stream, not just heat. |
