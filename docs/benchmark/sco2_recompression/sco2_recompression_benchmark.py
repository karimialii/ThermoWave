"""Benchmark: recompression supercritical-CO2 Brayton cycle vs. a real
published exergy analysis: M. Penkuhn and G. Tsatsaronis, "Exergoeconomic
analyses of different sCO2 cycle configurations", 6th Int. Symposium -
Supercritical CO2 Power Cycles, Pittsburgh, 2018. This is their
configuration (c), "Recompression, recuperated sCO2 cycle" (Figure 1c) --
NOT the fancier "modified recompression" (d), which has a third
compressor. All published numbers below are read from the paper's own
Table 5(c)/6(c).

This is the direct payoff of four capability gaps identified while
assessing a third-party TESPy reproduction of this same paper earlier
this session:
  #1 Junction free split fractions      -> the recompression splitter
  #2 Setpoint/Controller live targets   -> the T-match at the merge point
  #3 core/exergy.py                     -> the E_F/E_P/E_D comparison below
  #4 entropy on ideal-gas fluids        -> not used directly (CO2 uses
                                            CoolProp, which already had
                                            entropy), but SteamTurbine/
                                            Pump's entropy-based physics
                                            (needed for a REAL fluid near
                                            its critical point, where a
                                            gamma relation is badly wrong)
                                            is exactly the kind of thing
                                            those additions made routine.

CYCLE TOPOLOGY (matches sCO2.py's own wiring):
    Source(c4: 600 C, 250 bar)
      -> SteamTurbine (eta_s=0.90, P_out=77.95 bar)                    [c5]
      -> MultiPassHeatExchanger "rec2" hot side                       [c14]
      -> MultiPassHeatExchanger "rec1" hot side                       [c15]
      -> Junction "sp1" (free split -- THE recompression fraction)
           out0 -> Pipe "cooler" (fixed duty)          -> Pump "cp1"   [c2]
                -> MultiPassHeatExchanger "rec1" cold side            [c12]
           out1 -> Pump "cp2"                                         [c11]
      -> Junction "m1" (merge c12 + c11)                              [c13]
      -> MultiPassHeatExchanger "rec2" cold side                      [c3]
      -> Sink

Controller closes sp1's free split fraction by matching cp2's outlet
temperature to rec1's cold-outlet temperature -- exactly TESPy's own
`c11.set_attr(T=Ref(c12, 1, 0))` constraint, the real design rule that
mixing two streams at different temperatures at the merge destroys exergy
for no benefit.

WHY Pump AND SteamTurbine, NOT Compressor/SimpleCompressor/SimpleTurbine
--------------------------------------------------------------------------
SimpleCompressor/SimpleTurbine use the ideal-gas gamma relation
T2s = T1*PR**((gamma-1)/gamma) -- badly wrong this close to CO2's critical
point (31.1 C / 73.8 bar), where cp is many times its ideal-gas value (see
this session's own check: ~6000 J/kg-K near 35 C/75 bar vs ~850 J/kg-K
away from the dome). Pump and SteamTurbine are both entropy-based
(s_in -> h_out_isentropic via enthalpy_ps), fluid-model-agnostic despite
their liquid/steam-flavored names -- exactly the physics real-fluid
compression/expansion needs, regardless of what the working fluid is
called.

HOW THE BOUNDARY/DESIGN DATA WAS DERIVED
------------------------------------------
Component parameters (recuperator UA, mass flow, cooler duty) are NOT
invented -- they're derived once, outside ThermoWave, purely from the
paper's own Table 5(c) (T, p, e^T, e^M at every state) and Table 6(c)
(E_F/E_P/E_D/epsilon per component). Mass flow isn't published directly,
so it's back-derived from the exergy table -- but naively using the
turbine's published E_P (its exergy PRODUCT = shaft work) is wrong: SPECO
defines a turbine's E_P as shaft-side (ambiguous re: mechanical/generator
losses -- see Table 2's separate mechanical/generator efficiencies) while
its E_F is fluid-side (mdot*(e_in-e_out), mechanically unambiguous). Using
E_F for the turbine and E_P for each compressor (the fluid-side quantity
in both cases) gives three independent mass-flow estimates -- turbine,
cp1, cp2 -- that agree with each other to within 0.3%. See
sco2_recompression_benchmark.md for the full derivation and the resulting
comparison against the paper's own component efficiencies. The
recompression SPLIT ITSELF is genuinely solved by ThermoWave's own Newton
solve, not taken from the published data at all.

WHY TWO SOLVE STAGES
----------------------
Near-critical CO2's property surface is steep enough that a cold Newton
start with the split fraction ALSO free diverges into an unphysical state
(a finite-difference Jacobian step overshoots into a region CoolProp can't
even evaluate). Stage 1 solves the same network with the split fraction
FIXED at the derivation's own calibration value (no Controller needed --
see Junction's docstring: a single non-None entry alongside one None just
makes that one the fully-determined dependent fraction, not a free
unknown) to get every state close to right. Stage 2 then frees the split
and adds the Controller, warm-started from stage 1's converged result --
turning "solve the whole near-critical cycle cold" into "solve it, then
ask one more question of an already-good solution," which is exactly the
continuation strategy the solver's own convergence-failure message
recommends.

Run directly (needs the coolprop extra: pip install thermowave[coolprop]):
    python docs/benchmark/sco2_recompression/sco2_recompression_benchmark.py
"""

from __future__ import annotations

from thermowave.components import (
    Controller,
    Junction,
    MultiPassHeatExchanger,
    Pipe,
    Pump,
    Sensor,
    Sink,
    Source,
    SteamTurbine,
)
from thermowave.core import Network
from thermowave.core.exergy import exergy_report
from thermowave.fluids import CoolPropFluid

CO2 = CoolPropFluid(name="CO2")
BAR = 1.0e5

# --- Derived design point (see module docstring + .md for the derivation) --
# mdot_total/mdot_cp1/mdot_cp2 come from three independent, purely paper-
# based estimates (turbine E_F, cp1 E_P, cp2 E_P against Table 5(c)/6(c))
# that agree with each other to within 0.3% -- see the .md for the full
# derivation and why the earlier attempt (using the turbine's E_P instead
# of E_F) was off by ~2%.
MDOT_TOTAL = 1177.781         # kg/s
MDOT_CP1 = 857.314            # kg/s, cooler+cp1 branch
FRAC_CALIBRATION = MDOT_CP1 / MDOT_TOTAL   # 0.727906
Q_COOLER = 138.015e6          # W
UA_REC1 = 37_685.65e3         # W/K
UA_REC2 = 21_284.05e3         # W/K

# Pressures (bar), all read directly from connection_validation.csv
P4, P5 = 250.0, 77.95         # turbine in/out
P14, P15 = 76.94, 75.15       # rec1 hot out / rec2 hot out -> rec1 hot in/out
P1, P2 = 75.0, 258.4          # cp1 in/out
P10, P11 = 75.15, 257.51      # cp2 in/out
P12, P13, P3 = 257.51, 257.51, 257.0  # rec1 cold out / merge out / rec2 cold out

PR_HOT_REC2 = P14 / P5
PR_HOT_REC1 = P15 / P14
PR_COLD_REC1 = P12 / P2
PR_COLD_REC2 = P3 / P13
PR_CP1 = P2 / P1
PR_CP2 = P11 / P10


def build_network(split_fractions):
    """split_fractions=[float, float]: fully fixed (stage 1).
    split_fractions=[None, None]: sp1.frac0 free, closed by a Controller
    matching the merge temperatures (stage 2)."""
    net = Network(fluid=CO2)

    src = Source(name="c4", P=P4 * BAR, T=600.0 + 273.15, mdot=MDOT_TOTAL)
    turb = SteamTurbine(name="turb", P_out=P5 * BAR, eta_s=0.90)
    # counterflow effectiveness is exactly n_passes-invariant (see
    # MultiPassHeatExchanger's own docstring) -- n_passes=1 is the same
    # answer as any other value here, just far fewer internal nodes/faster.
    rec2 = MultiPassHeatExchanger(
        name="rec2", UA=UA_REC2, PR_hot=PR_HOT_REC2, PR_cold=PR_COLD_REC2,
        n_passes=1, arrangement="counterflow",
    )
    rec1 = MultiPassHeatExchanger(
        name="rec1", UA=UA_REC1, PR_hot=PR_HOT_REC1, PR_cold=PR_COLD_REC1,
        n_passes=1, arrangement="counterflow",
    )
    sp1 = Junction(name="sp1", n_inlets=1, n_outlets=2, split_fractions=split_fractions)
    cooler = Pipe(name="cooler", L=1.0, D=1.0, f=0.0, heat_loss=Q_COOLER)
    cp1 = Pump(name="cp1", PR=PR_CP1, eta=0.85)
    cp2 = Pump(name="cp2", PR=PR_CP2, eta=0.85)
    m1 = Junction(name="m1", n_inlets=2, n_outlets=1)
    sink = Sink(name="sink")

    components = dict(
        src=src, turb=turb, rec2=rec2, rec1=rec1, sp1=sp1, cooler=cooler,
        cp1=cp1, cp2=cp2, m1=m1, sink=sink,
    )

    is_free = bool(sp1.free_parameters())
    if is_free:
        sensor_cp2_out = Sensor(name="sensor_cp2_out")
        sensor_rec1_cold_out = Sensor(name="sensor_rec1_cold_out")
        match_recompression_split = Controller(
            name="match_recompression_split",
            sensor=sensor_cp2_out, quantity="T [K]",
            component=sp1, free_param="frac0",
            value=lambda s: sensor_rec1_cold_out.report_metrics(s)["T [K]"],
        )
        components.update(
            sensor_cp2_out=sensor_cp2_out,
            sensor_rec1_cold_out=sensor_rec1_cold_out,
            match_recompression_split=match_recompression_split,
        )

    for c in components.values():
        net.add_component(c)

    net.connect(src, "out", turb, "in")
    net.connect(turb, "out", rec2, "hot_in")
    net.connect(rec2, "hot_out", rec1, "hot_in")
    net.connect(rec1, "hot_out", sp1, "in0")
    net.connect(sp1, "out0", cooler, "in")
    net.connect(cooler, "out", cp1, "in")
    net.connect(cp1, "out", rec1, "cold_in")
    net.connect(sp1, "out1", cp2, "in")
    if is_free:
        net.connect(cp2, "out", components["sensor_cp2_out"], "tap")
        net.connect(components["sensor_cp2_out"], "tap", m1, "in1")
        net.connect(rec1, "cold_out", components["sensor_rec1_cold_out"], "tap")
        net.connect(components["sensor_rec1_cold_out"], "tap", m1, "in0")
    else:
        net.connect(cp2, "out", m1, "in1")
        net.connect(rec1, "cold_out", m1, "in0")
    net.connect(m1, "out0", rec2, "cold_in")
    net.connect(rec2, "cold_out", sink, "in")

    return net, components


print("--- stage 1: fixed split, establishing a good starting point ---")
net1, c1_components = build_network([FRAC_CALIBRATION, None])
result1 = net1.solve(tol=1e-3, max_iter=300, damping=0.5, step_growth=1.05, progress=False)
print(f"stage 1 converged: {result1.converged}, iterations: {result1.iterations}")

print("\n--- stage 2: free split, closed by the T-match Controller ---")
net2, c2 = build_network([None, None])
result2 = net2.solve(
    tol=1e-3, max_iter=300, damping=0.5, step_growth=1.05, progress=False, warm_start=result1
)
print(f"stage 2 converged: {result2.converged}, iterations: {result2.iterations}")

print("\n--- stage 3: polish to tight tolerance from stage 2's result ---")
result2 = net2.solve(
    tol=1e-8, max_iter=100, damping=0.6, step_growth=1.05, progress=False, warm_start=result2
)
print(f"stage 3 converged: {result2.converged}, iterations: {result2.iterations}")
result2.print_report()

state = result2.state()
frac0 = state.param("sp1.frac0")
print(f"\nsolved recompression split (frac to cp1 branch): {frac0:.4f}")
print(f"(derivation-calibration value used for UA/Q_cooler was {FRAC_CALIBRATION})")

T_cp2_out = c2["sensor_cp2_out"].report_metrics(state)["T [K]"] - 273.15
T_rec1_cold_out = c2["sensor_rec1_cold_out"].report_metrics(state)["T [K]"] - 273.15
print(f"T at merge: cp2 out = {T_cp2_out:.2f} C, rec1 cold out = {T_rec1_cold_out:.2f} C "
      f"(published: both 264.11 C)")

T3_solved = CO2.temperature_ph(*state.node("rec2.cold_out")) - 273.15
print(f"\npredicted rec2 cold outlet (c3): T = {T3_solved:.2f} C  (published: 433.63 C)")

turb, cp1, cp2 = c2["turb"], c2["cp1"], c2["cp2"]
report = exergy_report(
    result2, T0=298.15, P0=101325.0,
    fuel=[lambda s: MDOT_TOTAL * (CO2.enthalpy_pt(P4 * BAR, 600.0 + 273.15)
                                   - CO2.enthalpy_pt(P3 * BAR, 433.693 + 273.15))],
    product=[lambda s: turb.report_metrics(s)["power [W]"]
              - cp1.report_metrics(s)["power [W]"]
              - cp2.report_metrics(s)["power [W]"]],
)
report.print_report()

# --- comparison: ThermoWave vs. the paper's own Table 6(c). ---
# (component, E_F [MW], E_P [MW], E_D [MW], epsilon [%])
PUBLISHED = {
    "turb": (197.4, 185.3, 12.2, 93.8),
    "rec1": (73.9, 70.0, 3.9, 94.8),
    "rec2": (139.2, 135.4, 3.8, 97.3),
    "cp1": (47.6, 40.2, 7.3, 84.6),
    "cp2": (37.7, 32.9, 4.8, 87.3),
}
# Expected gap for turb/cp1/cp2, purely from the paper's own Table 2
# efficiency parameters: ThermoWave's Pump/SteamTurbine model only the
# working-fluid thermodynamics (no motor/generator), while the paper's
# Table 6(c) T-1/C-1A/C-1B rows appear to already fold those in.
EXPECTED_GAP_PTS = {
    "turb": (1 / 0.99 / 0.99 - 1) * 100,   # turbine mech. eff x generator eff
    "cp1": (1 / 0.98 / 0.97 - 1) * 100,    # compressor mech. eff x motor eff
    "cp2": (1 / 0.98 / 0.97 - 1) * 100,
    "rec1": None,  # mean-cp vs. inlet-cp UA convention gap, see the .md
    "rec2": None,
}

print(f"\n{'component':<10}{'TW eps':>9}{'pub eps':>9}{'delta':>9}{'expected':>10}")
print("-" * 47)
for name in ("turb", "rec1", "rec2", "cp1", "cp2"):
    cost = report.component_costs[name]
    eps = cost.epsilon * 100.0
    delta = eps - PUBLISHED[name][3]
    expected = EXPECTED_GAP_PTS[name]
    expected_str = f"~{expected:+.1f}" if expected is not None else "n/a"
    print(f"{name:<10}{eps:>9.2f}{PUBLISHED[name][3]:>9.2f}{delta:>+9.2f}{expected_str:>10}")

print("-" * 47)
print("\nSee sco2_recompression_benchmark.md for why turb/cp1/cp2's gaps match "
      "the paper's own Table 2 motor/mechanical/generator efficiencies, and "
      "why rec1/rec2's smaller gaps come from a UA-calibration convention "
      "difference instead.")
