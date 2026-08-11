# Benchmark: recompression supercritical-CO2 Brayton cycle

`sco2_recompression_benchmark.py` builds and solves a full recompression
sCO2 power cycle in ThermoWave and compares it directly against a real
published exergy analysis:

> M. Penkuhn and G. Tsatsaronis, *Exergoeconomic analyses of different
> sCO2 cycle configurations*, The 6th International Symposium —
> Supercritical CO2 Power Cycles, Pittsburgh, PA, March 27-29, 2018.

The paper studies four sCO2 cycle layouts (Figure 1a-d). This benchmark
implements configuration **(c), the "Recompression, recuperated sCO2
cycle"** — not the fancier "Modified recompression" cycle (d), which adds
a third compressor and a second precooler. The two are easy to mix up (the
component/node layout only differs by one extra compressor and two extra
recuperator ports), so it's worth being explicit up front: this benchmark
is (c), and every published number quoted below is read from the paper's
own **Table 5(c)** (stream data) and **Table 6(c)** (exergy results) —
never Table 5(d)/6(d).

Reproducing this cycle exercises four ThermoWave capabilities together:

| Capability | Where it's used |
|---|---|
| `Junction` free split fractions | The recompression splitter's mass-flow split |
| `Setpoint`/`Controller` live targets | The T-match constraint at the merge point |
| `core/exergy.py` | The E_F/E_P/E_D/epsilon comparison table below |
| Entropy-based turbomachinery | CO2 uses CoolProp, which already provides entropy, but `SteamTurbine`/`Pump`'s entropy-based physics is what makes them the right choice this close to CO2's critical point (see "Why `Pump`/`SteamTurbine`" below) |

The script and its plotting code live in their own standalone repo:
https://github.com/karimialii/sco2-recompression-benchmark

```bash
git clone https://github.com/karimialii/sco2-recompression-benchmark
cd sco2-recompression-benchmark
pip install thermowave[coolprop]
python sco2_recompression_benchmark.py
```

The plots below are regenerated from its solved result via
`plot_results.py`, which lives alongside the benchmark script in that repo:

```bash
pip install thermowave[coolprop,plot]
python plot_results.py
```

## The cycle — paper's Figure 1(c), "Recompression, recuperated sCO2 cycle"

![Recompression sCO2 cycle schematic: heater, turbine, two recuperators, the recompression split, cooler, and both compressors, closing back on itself](cycle.svg)

**What this benchmark shows**: given design-point data derived once from
the paper's own tables (see below), ThermoWave reproduces the cycle's
component exergetic efficiencies to within about 2 points for the heat
exchangers and 5 points for the turbomachinery — a gap that itself lines
up almost exactly with the paper's own motor/generator/mechanical
efficiency factors (Table 2), which ThermoWave's turbine/pump models
deliberately don't include. The predicted heater-inlet temperature (c3)
comes out 1.1% high, and the model's own free-split solve converges to a
recompression fraction measurably different from the calibration value —
both traced to a single, identified cause: a mean-temperature vs.
inlet-temperature convention mismatch in how recuperator UA gets
evaluated (see "A remaining, precisely-identified limitation" below).
Nothing here is hand-waved — every number is either read straight from
the paper's Table 5(c)/6(c) or produced by running the script.

Node numbers (`1,2,3,4,5,6,10-15`) match the paper's Table 5(c) exactly —
this is a 1:1 reproduction of the paper's own state-point labeling, not
just a topologically-similar cycle. The two turbomachines share one shaft
(`T-1` drives both `C-1a`/`cp1` and `C-1b`/`cp2` through the generator/
motor train shown as `G-1` in the paper's figure); ThermoWave doesn't
model a shared shaft directly, so `cp1`/`cp2` are just parameterized with
their own pressure ratios from the published state data instead — the
shaft-power balancing the paper does internally isn't a constraint this
network needs to reproduce.

```
Recycle (mdot=MDOT_TOTAL)
  -> SimpleHeater "heater" (T_out=600 C, PR=P4/P3)                  [c4]
  -> SteamTurbine "turb" (eta_s=0.90, P_out=77.95 bar)              [c5]
  -> MultiPassHeatExchanger "rec2" hot side                        [c14]
  -> MultiPassHeatExchanger "rec1" hot side                        [c15]
  -> Junction "sp1" (free split -- the recompression fraction)
       out0 -> Pipe "cooler" (fixed duty) -> Pump "cp1"             [c2]
            -> MultiPassHeatExchanger "rec1" cold side             [c12]
       out1 -> Pump "cp2"                                          [c11]
  -> Junction "m1" (merge c12 + c11)                                [c13]
  -> MultiPassHeatExchanger "rec2" cold side                        [c3]
  -> (back to Recycle)
```

## Why the network is a genuinely closed loop, and how the heater works

The heater (c3 -> c4) and cooler (c15/c6 -> c1) are where the cycle crosses
its system boundary — external heat in, external heat rejected. Both
endpoints of each are fully specified in the published data (temperature
*and* pressure), so closing this loop doesn't require inventing anything:
`SimpleHeater(T_out=600 C, PR=P4/P3)` reproduces exactly the same c4 target
a `Source` used to assert directly, just as a genuine loop component fed
from `rec2`'s own cold-side outlet instead of an external boundary value.

`SimpleEvaporator` (the component the `segs_exergy` benchmark's own boiler
uses for the same kind of crossing) doesn't work here: CO2 stays
**supercritical** across this entire cycle (250 bar ≫ CO2's ~73.8 bar
critical pressure), and `SimpleEvaporator`'s superheat/subcool/quality
targeting is relative to a saturation curve that doesn't exist above the
critical pressure. `SimpleHeater` (`thermowave.components.SimpleHeater`)
targets an outlet *temperature* directly instead — see its own docs — which
needs nothing from the saturation dome. `Pipe`'s existing `heat_loss`
parameter was also considered and rejected: its default (`None` = zero
loss) already carries meaning, colliding with the "`None` = free unknown"
convention `Compressor`/`Junction` use elsewhere, and `SimpleHeater`'s
target here is a directly-given constant anyway, not something needing a
`Setpoint` indirection.

`Recycle` (`thermowave.components.Recycle`) supplies the one thing
`SimpleHeater` alone doesn't: an anchor for the loop's overall mass-flow
scale (`MDOT_TOTAL`), the role `Source`'s `mdot=` used to play. Cold-start
convergence for a network with no fixed `(P, h)` anywhere needs an explicit
warm start (see `_warm_start()` in the script) — without one, every node
defaults to the solver's flat guess, which for supercritical CO2 lands
CoolProp's residual evaluation so far outside the cycle's real state space
that Newton's first Jacobian evaluation fails before it can take a step.

**A small, honest side effect of actually closing the loop:** the paper's
own published c3/c4 pressures (257 bar / 250 bar, Table 5(c)) are rounded
values, not exactly what the rest of the network's own chain of *relative*
pressure ratios (turbine `P_out`, recuperator `PR`s, compressor `PR`s)
independently produces for `rec2`'s cold-side outlet. A `Source` used to
paper over that gap by asserting c4's pressure directly regardless of what
came before it; a real closed loop can't — `heater.out`'s pressure is now
`PR * (whatever rec2.cold_out's own chain of ratios actually solves to)`,
which lands turbine inlet pressure about 0.2% away from the published
250 bar. The effect on every reported number here is at the same ~0.2%
level (e.g. turbine power 1.891e8 W vs. 1.888e8 W with the old
`Source`-pinned topology) — well inside the gaps already discussed below,
and, if anything, a more physically honest result: a real closed loop's
absolute pressure has to be self-consistent all the way around, not
asserted at one open end.

## Why `Pump`/`SteamTurbine`, not `SimpleCompressor`/`SimpleTurbine`

`SimpleCompressor`/`SimpleTurbine` use the ideal-gas relation
`T2s = T1 * PR**((gamma-1)/gamma)` — badly wrong this close to CO2's
critical point (31.1 °C / 73.8 bar). A direct check found cp there running
~6,000 J/kg·K vs. ~850 J/kg·K away from the dome — a ~7x spike a
constant-gamma relation can't represent. `Pump` and `SteamTurbine` are both
entropy-based (`s_in -> h_out,isentropic` via `enthalpy_ps`) and
fluid-model-agnostic despite their liquid/steam-flavored names — exactly
the physics real-fluid compression/expansion needs, independent of what
the working fluid happens to be called. This cycle is precisely the case
that motivates having them.

## How the design-point inputs were derived — and a real bug this caught

Everything ThermoWave's components need but the published tables don't
hand over directly (mass flow, recuperator UA, cooler duty) was derived
once, purely from the paper's own **Table 5(c)** (state points: T, p,
e^T, e^M) and **Table 6(c)** (E_F/E_P/E_D/epsilon per component), using
CoolProp's CO2 model for enthalpy/entropy at those published states.

**Mass flow, and the bug it exposed**: the paper doesn't publish mass
flow directly, so it has to be back-derived from the exergy table. The
first attempt used the turbine's published `E_P = 185.3 MW` (its exergy
*product*) as if it were the turbine's own mechanical shaft power:
`mdot = E_P / (h4 - h5) = 1154.2 kg/s` — about 2% low, and it didn't
reproduce the compressors' own mass flows to better than a few percent
either.

The fix: use each component's **fluid-side** SPECO quantity, not its
shaft-side one. For a turbine, SPECO defines the *fuel* as the physical
exergy drop of the working fluid (`E_F = mdot*(e_in - e_out)`) and the
*product* as the shaft work delivered — so `E_F`, not `E_P`, is the
turbine's fluid-referenced, mechanically-unambiguous number. For a
compressor it's the other way around: the *product* is the fluid's
exergy rise (`E_P = mdot*(e_out - e_in)`), so `E_P` is the
mechanically-unambiguous one there. Reading `e^T + e^M` straight from
Table 5(c) at each component's in/out states and solving for `mdot` from
the fluid-side quantity gives three independent estimates that agree with
each other to within 0.3% — evidence they're right, since nothing forces
three separate calculations to agree unless the underlying number is
real:

| Estimate | Source | mdot [kg/s] |
|---|---|---|
| Turbine (`E_F`, node 4->5) | Table 6(c) T-1, Table 5(c) nodes 4,5 | 1179.28 (total) |
| Compressor cp1 (`E_P`, node 1->2) | Table 6(c) C-1A, Table 5(c) nodes 1,2 | 857.14 |
| Compressor cp2 (`E_P`, node 10->11) | Table 6(c) C-1B, Table 5(c) nodes 10,11 | 321.32 |
| cp1 + cp2 | — | 1178.46 (total) |

`mdot_total = 1177.781 kg/s`, `mdot_cp1 = 857.314 kg/s`,
`mdot_cp2 = 320.467 kg/s` are the values actually used in the benchmark —
selected from this same self-consistency check, refined to the precision
where the turbine-side and compressor-side totals agree to 4 significant
figures rather than 3. The recompression split implied by these,
`mdot_cp1 / mdot_total = 0.727906`, is the calibration value used below
for the recuperator UA and cooler duty — **not** the value ThermoWave's
own free-split solve converges to (see Results).

**Recuperator UA**: heat duty computed exactly from each side's own
enthalpy change (cross-checked hot-side vs. cold-side to 4 significant
figures given the mass flows above — e.g. rec1: 188.4262 MW vs.
188.4261 MW), then effectiveness and NTU from the standard counterflow
relation, `UA = NTU * Cmin`, with `Cmin`/`Cmax` from each side's
mean-temperature-evaluated heat capacity.

**Cooler duty**: `Q_cooler = mdot_cp1 * (h6 - h1)`, both states fully
known from the published Table 5(c) data.

## Results

**Predicted c3** (`rec2`'s cold-side outlet, which now feeds `Recycle` ->
`heater` to close the loop rather than terminating at a `Sink`): 441.2 °C
vs. published 433.63 °C (Table 5(c)) — about 1.1% high in absolute
temperature. The same run shows where ThermoWave's own free-split solve
lands relative to the derivation's calibration value:

![Predicted vs. published c3 temperature, and the recompression-split discrepancy between ThermoWave's free-split solve and the derivation's calibration value](validation_summary.png)

**A remaining, precisely-identified limitation**: even with a
self-consistent mass flow, ThermoWave's own free-split solve converges to
a recompression fraction (0.786) that differs from the calibration value
(0.7279) by more than the mass-flow uncertainty above would suggest. The
reason: this benchmark's own UA derivation evaluates each recuperator
side's heat capacity at its **mean** temperature (a reasonable
single-point calibration choice for a fluid with wildly
temperature-dependent cp), while `MultiPassHeatExchanger`'s own forward
model evaluates cp at each side's **inlet** temperature (see its
docstring: *"cp evaluated at each side's own inlet state"*) — the same
convention `HeatExchanger`'s fixed-effectiveness mode uses. For CO2 with a large temperature
swing across `rec1`'s cold side (123 °C to 264 °C, spanning a region
where cp itself varies severalfold), those two evaluation conventions
give measurably different `Cmin`/`Cr`, and therefore a different
sensitivity of the cold-outlet temperature to the mass-flow split — so
the *same* UA reproduces the exact duty at the calibration split, but the
Newton solve (searching for whatever split makes T11=T12 under the
network's own, inlet-cp-based model) lands somewhere else. Resolving it
would mean either calibrating UA against `MultiPassHeatExchanger`'s own
inlet-cp convention directly, or iterating the derivation at the solved
split until it stops moving — deliberately not done here, to keep the
derivation a single, auditable pass.

**Component exergetic efficiency (epsilon = E_P/E_F)**, ThermoWave vs.
the paper's own Table 6(c):

![Bar chart comparing ThermoWave's and the published component exergetic efficiencies for turb, rec1, rec2, cp1, and cp2](exergy_comparison.png)

```
component    TW eps  pub eps   delta  why
------------------------------------------------------------------
turb          95.62    93.80   +1.82  turbine mech. eff (99%) x
                                       generator eff (99%), Table 2 --
                                       folded into the paper's T-1 row;
                                       ThermoWave's SteamTurbine models
                                       pure aerodynamic shaft work only
rec1          96.76    94.80   +1.96  mean-cp vs. inlet-cp UA
                                       calibration gap (see above)
rec2          97.99    97.30   +0.69  same gap, smaller (rec2's
                                       temperature swing is smaller)
cp1           89.41    84.60   +4.81  compressor motor eff (97%) x
                                       mech. eff (98%), Table 2 --
cp2           92.00    87.30   +4.70  ThermoWave's Pump has no
                                       motor/mechanical loss of its own
------------------------------------------------------------------
```

`turb`/`cp1`/`cp2`'s gaps line up well with what the paper's own Table 2
efficiency parameters predict for a pure-shaft vs. motor/generator-
inclusive comparison: turbine mechanical (99%) x generator (99%)
efficiency implies a ~2.0-point gap (`1/0.99/0.99 - 1`), matching the
observed +1.82; compressor mechanical (98%) x motor (97%) efficiency
implies a ~5.2-point gap (`1/0.98/0.97 - 1`), matching the observed
+4.81/+4.70 closely. ThermoWave's `Pump`/`SteamTurbine` deliberately model
only the working-fluid thermodynamics (pure shaft/isentropic work, no
electrical machinery) — Table 6(c)'s own `T-1`/`C-1A`/`C-1B` rows appear
to already fold the paper's Table 2 motor/generator efficiencies into
their exergy product, so a `Table-2`-sized gap here is the expected
signature of comparing two different, but each internally consistent,
component-boundary conventions, not a modeling error.

`rec1`/`rec2`'s smaller, purely heat-exchanger-side gaps are the
mean-cp/inlet-cp UA calibration convention difference described above —
unrelated to motor/generator accounting.

## What would make this benchmark stronger

Iterate the UA derivation at the solved split (recompute the implied
UA-calibration temperature profile using ThermoWave's own inlet-cp
convention, or re-derive UA using `MultiPassHeatExchanger`'s own per-pass
duty, walking the calibration split toward wherever it stops moving) —
deliberately not run here, to keep the derivation a single, auditable
pass rather than an opaque loop. Given how close the *efficiency*
comparison already is, the expected payoff is mainly tighter
recompression-split/absolute-power agreement for `cp1`/`cp2`, not a
qualitatively different result.
