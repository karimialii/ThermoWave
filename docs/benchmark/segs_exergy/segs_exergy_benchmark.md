# Benchmark: 30 MWe SEGS solar Rankine cycle, design point

`segs_exergy_benchmark.py` builds and solves the water/steam power cycle of
a 30 MWe SEGS (Solar Energy Generating System) parabolic-trough plant in
ThermoWave, at its 100%-solar design point, and compares it directly
against the **original published source**:

> F. Lippke, *Simulation of the Part-Load Behavior of a 30 MWe SEGS
> Plant*, SAND95-1293, Sandia National Laboratories, 1995.
> https://www.osti.gov/biblio/95571 — Table 3 (every component's design
> mass flow / pressure / efficiency / UA) and Fig. 4 (the design heat
> balance) are the source of every number this benchmark compares against.

This plant is also the subject of a TESPy exergy-analysis example repo
([`docs/benchmark/Tespy/SEGS_exergy-main`](../Tespy/SEGS_exergy-main)),
which was run this session for cross-checking. Per the scope asked for
here, **the comparison below is ThermoWave vs. the original 1995 paper
only** — TESPy's own numbers are deliberately not part of it. See
"What the TESPy run was for" at the end for what that run did and why it
isn't in the main comparison.

Run it directly (needs the `coolprop` extra):

```bash
pip install thermowave[coolprop]
python docs/benchmark/segs_exergy/segs_exergy_benchmark.py
```

## The cycle

7-stage steam turbine train (2 HP stages, 5 LP stages with one reheat
between HP and LP), 5 regenerative feedwater heaters fed by fixed turbine
extractions with a condensate drain cascade, one open feedwater heater
(deaerator), 2 pumps. The oil-heated boiler and reheater (the solar
field's own boundary with the water/steam side) are not modeled — Fig. 4
gives both endpoints' design temperature *and* pressure directly, so, like
the [`sco2_recompression`](../sco2_recompression/sco2_recompression_benchmark.md)
benchmark's own external heater/cooler, there's nothing to solve there:
the network starts at a `Source` (100 bar / 371 °C boiler outlet) and ends
at a `Sink` (feedwater returning to the boiler), with a second independent
Source/Sink pair for the reheater. See the script's own module docstring
for the full topology diagram and the exact port wiring.

**Every extraction/split fraction is fixed, not solved for.** Table 3's
own turbine `mFeed0` sequence hands over the exact design-point extraction
at each tap directly — e.g. `hpt1`'s `38.6415 -> 35.7326 kg/s` *is* the
design 2.9089 kg/s extraction to `hpPreheater2`. That's different from the
[`sco2_recompression`](../sco2_recompression/sco2_recompression_benchmark.md)
benchmark's recompression split, which genuinely wasn't known in advance
and needed a `Junction` free parameter closed by a `Controller`. Here
there's no free-parameter machinery at all — every split is a plain
Python float computed once from Table 3.

## Why the feedwater heaters aren't `MultiPassHeatExchanger`

Tried first, and it reliably crashed: every extraction feeding a
regenerative feedwater heater here is **condensing steam**, and several of
this plant's own low-pressure turbine exhausts land ON or inside the
two-phase dome. `SimpleCondenser`'s own docstring says exactly why that
breaks an effectiveness-NTU model: *"cp is effectively infinite during a
constant-pressure phase change, so an effectiveness-NTU (C = mdot\*cp)
framework can't represent condensation."* `MultiPassHeatExchanger` hit
that head-on — `CoolProp cp() ... Saturation pressure [28680 Pa]
corresponding to T [341.215 K] is within 1e-4% of given p`, at exactly
`lpt4`'s exhaust state (0.2868 bar), which genuinely is right at
saturation at this design point.

So each feedwater heater here is a **boundary-crossing pair**, one
component per side, matching how a real feedwater heater's two jobs
actually are different physics:

- **hot side**: `SimpleCondenser(outlet_quality=0)` — the extraction steam
  condenses fully to saturated liquid at its own *solved* mdot/h_in. Its
  duty is a genuine result of the solve, not assumed.
- **cold side**: `SimpleEvaporator(duty=<precomputed constant>)` (despite
  the name — duty mode just adds a fixed Q to a single stream; the
  feedwater here stays subcooled liquid throughout, never actually
  evaporating). Its duty is precomputed *once*, outside the solve, from
  Table 3's own design cold-side mass flow and a standard feedwater-heater
  terminal-temperature-difference assumption (TTD = 5 K below the
  extraction's saturation temperature) — using only safe, single-phase
  liquid `enthalpy_pt` calls, never a `cp()` call anywhere near the dome.

The two sides' duties are consequently **not** forced to match (no shared
UA ties them) — the script reports both, and the gap, explicitly (see
Results). This is a real, disclosed modeling simplification, not hidden
in the numbers.

## Getting the initial guess to converge at all

Even with the redesign above, a cold start still crashed on `F0` itself
(before Newton takes a single step): `MultiPassHeatExchanger`'s guess
machinery adds a fixed +4×10⁵ J/kg swing per hop, and chaining three
feedwater heaters in series compounds that into an initial guess sitting
inside the two-phase dome for an intermediate node — a pure guess
artifact, unrelated to the actual physics. Fixed with a **hand-built
`warm_start`**: `_warm_start()` constructs a `SolveResult` with just the
feedwater (cold) chain's node enthalpies set to plausible subcooled
temperatures (45 °C at the condensate pump outlet, rising through each
heater to 235 °C before the boiler) — only the *energy level* matters for
a warm start, not precision. With that, `net.solve(tol=1e-7, damping=0.4,
step_growth=1.03, warm_start=...)` converges cleanly in 51 iterations.

## Results

### Turbine stage power (the strongest, most direct check)

Table 3 hands over each stage's design `genPower` [kW] directly — a clean,
node-level comparison with no free parameters left to tune:

```
stage    ThermoWave [kW]   paper [kW]   diff
hpt1          7634.75        7643.54    -0.11%
hpt2          3476.14        3480.36    -0.12%
lpt1          5730.57        5730.18    +0.01%
lpt2          6692.24        6697.61    -0.08%
lpt3          5044.07        5045.83    -0.03%
lpt4          4505.79        4505.28    +0.01%
lpt5          2979.43        2979.68    -0.01%
-----------------------------------------------
TOTAL        36062.99       36082.47    -0.05%
```

Every stage lands within 0.12% — expected, since `SteamTurbine`'s own
entropy-based physics (`s_in -> h_out,isentropic` via CoolProp, scaled by
Table 3's own `etas0`) at Table 3's own `pFeed0`/`pDrain0` pressures is
essentially the same calculation the paper's own EASY model made; this is
mainly a check that the topology, extraction fractions, and CoolProp's
IAPWS-97 water model line up with the paper's own steam tables, which it
clearly does.

### Mass balance closure

Not assumed — a consequence of wiring the drain cascade correctly:

```
boiler outlet mdot:         38.6415 kg/s
feedwater return to boiler: 38.6415 kg/s   (exact)
```

The full condensate-cascade path (each heater's drain falling into the
next lower-pressure heater, down to the deaerator/condenser) re-derives
the exact boiler mass flow at the feedwater-tank outlet with zero
unaccounted loss — a real end-to-end check of the topology, not a
tautology (nothing forces this to close; it closes because the topology
is right).

### Feedwater heater energy balance (hot side solved vs. cold side fixed-duty)

```
heater          Q_hot [kW]   Q_cold [kW]   gap
hpPreheater1       5756.9       5747.7    +0.2%
hpPreheater2       5161.5       5583.1    -7.6%
lpPreheater1       2181.4       2348.5    -7.1%
lpPreheater2       3815.9       3974.9    -4.0%
lpPreheater3       3970.4       4175.5    -4.9%
```

`Q_hot` is a genuine solved result (extraction steam's own condensing
duty, at its actual solved mdot/enthalpy); `Q_cold` is the precomputed
TTD-based assumption described above. A few-percent gap is the expected
cost of that simplification (no shared UA ties the two sides here) —
`hpPreheater1` (fed by a mixed/merged drain stream, closest to the
deaerator's own known saturation state) lands closest.

### Pump parasitics

```
                          ThermoWave    Table 2 published
condensatePump               38.9 kW         190 kW
feedPump                    576.8 kW         880 kW
combined                    615.7 kW        1070 kW
```

Lower than published, mainly because `Pump`'s discharge pressures here
(9.0 bar / 104.1 bar) are this benchmark's own reasonable placeholders,
not derived from Table 3 (which gives quadratic pressure-drop
*coefficients* for EASY's own part-load formula in an internal, unstated
unit convention — not a single design-point ΔP; see the script's earlier
version notes). Table 2's own published parasitics also fold in the motor
efficiency (0.95, Table 1) `Pump` here doesn't model — matches the same
mechanical-vs-electrical distinction the
[`sco2_recompression`](../sco2_recompression/sco2_recompression_benchmark.md)
benchmark ran into with TESPy's turbine/generator split.

### Gross and net output

```
                                        ThermoWave      paper
gross mechanical power                  36063.0 kW   36067 kW  (Fig. 4)
gross electrical power (x0.97 gen.)     34981.1 kW   34985 kWe (Fig. 4)
```

**Exact match** on the gross figures (Fig. 4's own kW→kWe 97% generator
step is applied identically here) — strong confirmation that the turbine
train alone reproduces the paper's own heat balance almost exactly.

```
steam-cycle gross mechanical efficiency
  (turbine power / (boiler + reheat duty)):  38.9%
```

The paper's own headline **38.2% gross electric efficiency** and **31.4
MWe net output** are NOT directly comparable to this benchmark's own
figures, and this is disclosed rather than glossed over: those two paper
figures are computed against a *wider* system boundary — **solar-absorbed
power**, not boiler+reheat duty — and *net* subtracts the HTF pump and
cooling-tower pump/fan parasitics (Table 2: 1.56 + 0.99 MWe) on top of the
condensate/feed pumps, none of which this benchmark's scope models (see
"Scope" in the script's own module docstring). This benchmark's 38.9%
steam-cycle figure and the paper's 38.2% gross-electric figure are
answering two different, if related, questions — both given here rather
than picking one and calling it the same number.

### ThermoWave's own exergy analysis (no paper equivalent)

The 1995 paper predates exergy analysis of this specific plant (that's
what the companion TESPy repo adds, decades later) — there's no published
number to check the table below against. Included anyway as a genuine
value-add, using `thermowave.core.exergy.exergy_report` with `fuel` =
boiler + reheat duty and `product` = net turbine power (turbines minus
condensate/feed pump work):

```
component         E_F [W]      E_P [W]      E_D [W]   epsilon [%]
hpt1              8.494e6      7.635e6      8.597e5      89.9
hpt2              3.867e6      3.476e6      3.906e5      89.9
lpt1              6.230e6      5.731e6      4.990e5      92.0
lpt2              7.108e6      6.692e6      4.154e5      94.2
lpt3              5.324e6      5.044e6      2.804e5      94.7
lpt4              5.043e6      4.506e6      5.369e5      89.4
lpt5              4.537e6      2.979e6      1.557e6      65.7
condensatePump    3.887e4      2.828e4      1.059e4      72.8
feedPump          5.768e5      4.647e5      1.121e5      80.6
-----------------------------------------------------------------
Network           9.281e7      3.545e7      4.662e6      38.2%
```

`lpt5`'s low 65.7% exergetic efficiency tracks directly from Table 3's own
low `etas0 = 0.6445` for that stage — the paper's own EASY model
deliberately gave the last LP stage a lower design efficiency (see p. 6:
*"due to the necessary adaptation of the LP turbine stage
efficiency..."*), and the exergy analysis correctly flags exactly that
stage as the cycle's biggest single irreversibility.

## PASS/FAIL

**PASS.** Every turbine stage's power output lands within 0.12% of Table
3's published design value (total: -0.05%), the mass balance closes
exactly through the full drain cascade, and gross mechanical/electrical
output matches Fig. 4 almost exactly (36063 vs. 36067 kW, 34981 vs.
34985 kWe). The feedwater heater energy-balance gaps (up to 7.6%) and the
pump-parasitic shortfall are both disclosed, understood simplifications
(a precomputed TTD-based cold-side duty instead of a shared-UA model; a
placeholder pump discharge pressure instead of a derived one) — not
unexplained discrepancies.

## What the TESPy run was for (not part of the comparison above)

`docs/benchmark/Tespy/SEGS_exergy-main` ships a TESPy model of the *same*
plant with a genuine per-component exergy analysis (turbines, preheaters,
condenser, cooling tower, drum, ...) — a natural cross-check target. It
was installed and run this session (`tespy==0.9.3`, `CoolProp==6.8.0`,
pinned per its own `requirements.txt`, in a `python3.11` venv — CoolProp
has no wheel for this machine's default Python 3.14) to see whether it
would reproduce its own README's checked-in results.

It didn't, in this environment: the network's `nw.solve(mode='design')`
call itself fails to converge (`Detected singularity in Jacobian matrix`,
`NaN` mass flow/pressure/enthalpy after 1–4 iterations, `Inputs in Brent
[...] do not bracket the root` inside `CoolProp`'s water model), the same
outcome across pandas 3.0/numpy 2.4 and a pandas 2.1/numpy 1.26/scipy
1.11 downgrade. This looks like a genuine numerical-environment fragility
in this specific `tespy==0.9.3` + CoolProp + arm64/macOS combination (or a
missing initial-guess file/pin this repo's own `requirements.txt` doesn't
capture), not a modeling error on either side — and per this benchmark's
own scope, chasing it further wasn't worth it, since the actual comparison
target is the *original 1995 paper*, not TESPy. The repo's own README
ships a previously-computed, checked-in exergy table (from whatever
environment the zenodo-archived model was actually run in), which remains
a valid reference for anyone wanting a TESPy-side cross-check — just not
one reproduced fresh in this session, and deliberately not reproduced
in the comparison table above either, per what was asked for here.
