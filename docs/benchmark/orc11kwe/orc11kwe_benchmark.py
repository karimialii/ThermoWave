"""Benchmark: 11 kWe ORC waste-heat-recovery rig vs. ThermoWave's
Source -> Pump -> SimpleEvaporator chain.

An 11 kWe organic Rankine cycle (ORC) rig (R245fa working fluid, plate heat
exchangers, screw expander at a fixed 5000 rpm) was built and instrumented
under the ORCNext project to produce a public steady-state validation
dataset. The paper reports the *achieved operating envelope* across a
two-level design of experiments -- two heat-source mass flow rates
(1.5 kg/s and 3.0 kg/s) crossed with two heat-source inlet temperatures
(110 C and 120 C) -- as a min/max range table (their Table 5), not a
single design point:

    Variable                  Min. value     Max. value
    mdot_hf        [kg/s]     1.495          3.006
    T_hf,in        [C]        110.0          120.0
    T_hf,out       [C]        83.0           103.4
    p_evap         [bar]      9.510          12.457
    p_cond         [bar]      2.148          2.771
    T_superheat    [C]        19.8           21.1
    mdot_wf        [kg/s]     0.2902         0.3874

  Source: "Organic Rankine Cycle Part-Load Characterization: Validated
  Models of an 11 kWe ORC", University of Pretoria repository (ORCNext
  project dataset), Table 5.

  Note on p_evap/p_cond: pdftotext's OCR of the source PDF mangled both
  pressure variables' subscripts to the same "evap,in"/"evap,out" text
  (a rendering artifact of a Unicode-subscript font, not a paper error).
  The magnitudes make the physical identity unambiguous: R245fa's
  saturation pressure at ~100-120 C (the evaporator side) is in the
  9-12 bar range, and at ~25-35 C (the cooling-water condenser side) is in
  the 2-3 bar range -- so the second column is read here as the condensing
  pressure p_cond, not a second evaporator-side pressure.

WHAT THIS BENCHMARK DOES AND DOES NOT CHECK
---------------------------------------------
The paper validates net electrical power to within +/-2% and reports
evaporator/condenser heat-balance closure (primary-side duty vs.
secondary-side duty) within +/-5% -- but the specific numeric power/duty
values live in scatter plots (their Figures 5-9), not in a table, and the
paper's own expander/pump isentropic efficiencies aren't in the excerpt
available here either. Reproducing the +/-2% net-power figure exactly
would mean inventing an expander efficiency to hit a target, which is
fitting, not validation -- so this script does not attempt it.

What IS fully grounded in Table 5 (no invented numbers) is the *evaporator
heat balance*: the hot-fluid (Therminol 66 thermal oil) duty implied by
its published mdot_hf/T_hf,in/T_hf,out, versus the working-fluid duty
ThermoWave computes for the same evaporator (mdot_wf, p_evap, superheat).
Comparing those two against the paper's own claimed +/-5% closure band is
a real, citable check, not a tuned fit -- exactly the case the paper's own
Figures 8-9 make for their model.

The oil-side duty needs oil cp(T), which isn't published in the paper
excerpt either; it comes from Eastman's own published Therminol 66
correlation instead (see THERMINOL_66_CP below) -- a legitimate, citable
number, just from a different public source than the ORC paper itself.

The only genuinely assumed (not published, not derived) number is the
pump's isentropic efficiency, PUMP_ETA_ASSUMED. It matters very little
here: pump work for a liquid is dh ~= v*dP/eta, and for R245fa across this
cycle's pressure rise that is a few hundred J/kg against an evaporator duty
of ~2x10^5 J/kg -- under 1% of the quantity actually being checked. It is
included anyway (rather than skipped) so the network models a physically
complete pump-to-evaporator chain instead of silently pretending the pump
has zero work, and is called out explicitly here and in
orc11kwe_benchmark.md so nobody mistakes it for a validated figure.

Run directly (needs the coolprop extra: pip install thermowave[coolprop]):
    python docs/benchmark/orc11kwe/orc11kwe_benchmark.py
"""

from __future__ import annotations

from thermowave.components import Pump, Sink, SimpleEvaporator, Source
from thermowave.core import Network
from thermowave.fluids import CoolPropFluid

R245FA = CoolPropFluid(name="R245fa")

BAR = 1.0e5  # Pa

# --- Published operating envelope (Table 5) ---------------------------------
# The two columns are read as the two extreme corners of the paper's own
# 2x2 design of experiments (low mdot_hf/low T_hf,in vs. high/high) -- a
# reasonable pairing given how the DOE is described in the text, though the
# paper's Table 5 itself only states them as the achieved min/max range
# across all runs, not as two individually-paired rows.
CASES = {
    "low (mdot_hf~1.5 kg/s, T_hf,in=110 C)": dict(
        mdot_hf=1.495, T_hf_in_C=110.0, T_hf_out_C=83.0,
        p_evap=9.510 * BAR, p_cond=2.148 * BAR,
        T_superheat=19.8, mdot_wf=0.2902,
    ),
    "high (mdot_hf~3.0 kg/s, T_hf,in=120 C)": dict(
        mdot_hf=3.006, T_hf_in_C=120.0, T_hf_out_C=103.4,
        p_evap=12.457 * BAR, p_cond=2.771 * BAR,
        T_superheat=21.1, mdot_wf=0.3874,
    ),
}

# Eastman Therminol 66 published liquid-heat-capacity correlation
# (Cp in kJ/kg-K, T in deg C); Eastman TF-8695 technical bulletin.
def therminol66_cp(T_C: float) -> float:
    return (1.496005 + 0.003313 * T_C + 0.0000008970785 * T_C**2) * 1000.0  # J/kg-K

PUMP_ETA_ASSUMED = 0.75  # NOT published -- see module docstring
HEAT_BALANCE_TOL = 0.05  # the paper's own reported closure band

print(f"{'case':<40}{'Q_wf [kW]':>12}{'Q_oil [kW]':>12}{'rel. dev':>10}   result")
print("-" * 90)

all_within_tol = True
for label, c in CASES.items():
    net = Network(fluid=R245FA)

    T_sat_cond = R245FA.saturation_temperature(c["p_cond"])
    # 0.5 K of subcooling is a numerical-safety margin (keeps the source
    # state unambiguously single-phase liquid for CoolProp), not a claim
    # about the rig's actual subcooling.
    src = Source(name="src", P=c["p_cond"], T=T_sat_cond - 0.5, mdot=c["mdot_wf"])
    pump = Pump(name="pump", P_out=c["p_evap"], eta=PUMP_ETA_ASSUMED)
    evap = SimpleEvaporator(name="evap", superheat=c["T_superheat"])
    sink = Sink(name="sink")

    for comp in (src, pump, evap, sink):
        net.add_component(comp)
    net.connect(src, "out", pump, "in")
    net.connect(pump, "out", evap, "in")
    net.connect(evap, "out", sink, "in")

    result = net.solve(progress=False)
    assert result.converged, f"{label}: network did not converge"

    h_src = result.node_h["src.out"]
    h_evap_out = result.node_h["evap.out"]
    Q_wf = c["mdot_wf"] * (h_evap_out - h_src)  # W, working-fluid-side duty

    T_mean_oil_C = (c["T_hf_in_C"] + c["T_hf_out_C"]) / 2.0
    cp_oil = therminol66_cp(T_mean_oil_C)
    Q_oil = c["mdot_hf"] * cp_oil * (c["T_hf_in_C"] - c["T_hf_out_C"])  # W, oil-side duty

    rel_dev = (Q_wf - Q_oil) / Q_oil
    ok = abs(rel_dev) < HEAT_BALANCE_TOL
    all_within_tol &= ok
    print(f"{label:<40}{Q_wf / 1e3:>12.2f}{Q_oil / 1e3:>12.2f}{rel_dev:>+10.1%}   "
          f"{'PASS' if ok else 'FAIL'}")

print("-" * 90)
print(f"\n{'PASS' if all_within_tol else 'FAIL'}: ThermoWave's evaporator-side duty "
      f"(Source -> Pump -> SimpleEvaporator) matches the oil-side duty implied by the "
      f"published mdot_hf/T_hf,in/T_hf,out to within the paper's own reported "
      f"+/-{HEAT_BALANCE_TOL:.0%} heat-balance closure band, at both ends of the "
      f"published operating envelope.")

assert all_within_tol, "evaporator heat balance fell outside the paper's own +/-5% band"
