"""Benchmark: NASA Rotor 37 design point vs. ThermoWave's SimpleCompressor.

NASA Rotor 37 is a transonic axial compressor rotor designed and tested by
Reid and Moore at NASA Lewis (now Glenn Research Center), and has been used
as an open CFD/turbomachinery validation case since the 1994 IGTI blind test.
Its design-point aero-thermodynamics are public:

    Design mass flow            20.19 kg/s
    Total pressure ratio        2.106
    Adiabatic efficiency        0.877 (87.7%)
    Rotational speed            17,188.7 rpm
    Rotor tip speed             454.14 m/s
    Inlet total pressure        101,325 Pa  (standard day)
    Inlet total temperature     288.15 K    (standard day)

  Source: Ameri, "NASA Rotor 37 CFD Code Validation - Glenn-HT Code",
  NASA/CR-2010-216235 (AIAA-2009-1060), Table I, design parameters
  reproduced from Suder, NASA TM 107310 (1994).

WHAT THIS BENCHMARK DOES AND DOES NOT CHECK
---------------------------------------------
NASA's own validation of Rotor 37 is 3D CFD against blade-resolved flow-field
measurements (spanwise pressure/temperature profiles, tip vortex structure,
etc.) -- that is not what a 1D network solver can reproduce, and ThermoWave
does not attempt to.

What ThermoWave *can* validate is the aero-thermodynamic bookkeeping: given
the rotor's published (PR, eta_s) at its design mass flow and inlet state,
does the isentropic-compression-with-efficiency relation that
SimpleCompressor implements reproduce the exit temperature and shaft power
implied by those same published numbers? This script computes that expected
exit state and power by hand from the textbook relations, then checks
ThermoWave's own solved network against that hand calculation.

This is a self-consistency check on the thermodynamic relations (residual
correctness), anchored to NASA's published operating point, not a
reproduction of the internal flow field or the flow field non-uniformity
NASA measured. See rotor37_benchmark.md for the full discussion, including
why a proper multi-speed-line map validation (via
thermowave.components.Compressor) is future work, not done here: NASA/CR-
2010-216235 only tabulates the choke-to-stall sweep as figures, not a
digitized table, so there is no published multi-point map data to build a
.cop file from without manually reading points off a plot.

Run directly: `python docs/benchmark/rotor37/rotor37_benchmark.py`
(no optional extras needed -- IdealGasFluid has no external dependency).
"""

from __future__ import annotations

import math

from thermowave.components import SimpleCompressor, Sink, Source
from thermowave.core import Network
from thermowave.fluids import IdealGasFluid

# --- Published NASA Rotor 37 design point -----------------------------------
# Ameri, NASA/CR-2010-216235, Table I (design parameters, reproduced from
# Suder, NASA TM 107310).
MDOT_DESIGN = 20.19      # kg/s
PR_DESIGN = 2.106        # total-to-total pressure ratio
ETA_DESIGN = 0.877       # adiabatic (isentropic) efficiency
N_DESIGN_RPM = 17_188.7  # rev/min -- reported for context, not used below:
                          # SimpleCompressor has no map/shaft-speed coupling,
                          # see the module docstring.

# Standard-day inlet total conditions used for the NASA Lewis/Glenn rig.
P_IN = 101_325.0   # Pa
T_IN = 288.15      # K

GAMMA_AIR = 1.4  # standard value for the "textbook ideal-gas" comparison

# --- Fluid + network ---------------------------------------------------------
air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
net = Network(fluid=air)

src = Source(name="src", P=P_IN, T=T_IN, mdot=MDOT_DESIGN)
rotor37 = SimpleCompressor(name="rotor37", PR=PR_DESIGN, eta_s=ETA_DESIGN, gamma=GAMMA_AIR)
sink = Sink(name="sink")

net.add_component(src)
net.add_component(rotor37)
net.add_component(sink)
net.connect(src, "out", rotor37, "in")
net.connect(rotor37, "out", sink, "in")

result = net.solve()
result.print_report()

metrics = rotor37.report_metrics(result.state())
P_out, h_out = result.node_P[f"{rotor37.name}.out"], result.node_h[f"{rotor37.name}.out"]
T_out_solved = air.temperature_ph(P_out, h_out)
power_solved = metrics["power [W]"]

# --- Hand calculation, from the same published (PR, eta_s) -----------------
# T2s = T1 * PR ** ((gamma-1)/gamma)               (isentropic relation)
# dT_actual = (T2s - T1) / eta_s                    (efficiency definition)
# power = mdot * cp * dT_actual
T_out_isentropic = T_IN * PR_DESIGN ** ((GAMMA_AIR - 1.0) / GAMMA_AIR)
dT_isentropic = T_out_isentropic - T_IN
dT_actual = dT_isentropic / ETA_DESIGN
T_out_hand = T_IN + dT_actual
power_hand = MDOT_DESIGN * 1005.0 * dT_actual

print("\n--- Rotor 37 design-point benchmark -----------------------------")
print(f"{'quantity':<28}{'ThermoWave':>15}{'hand calc':>15}{'units':>8}")
print(f"{'PR':<28}{metrics['PR [-]']:>15.4f}{PR_DESIGN:>15.4f}{'-':>8}")
print(f"{'eta_s':<28}{metrics['eta_s [-]']:>15.4f}{ETA_DESIGN:>15.4f}{'-':>8}")
print(f"{'T_out':<28}{T_out_solved:>15.3f}{T_out_hand:>15.3f}{'K':>8}")
print(f"{'power':<28}{power_solved / 1e6:>15.4f}{power_hand / 1e6:>15.4f}{'MW':>8}")

tol = 1e-6
checks = [
    ("PR", metrics["PR [-]"], PR_DESIGN),
    ("eta_s", metrics["eta_s [-]"], ETA_DESIGN),
    ("T_out", T_out_solved, T_out_hand),
    ("power", power_solved, power_hand),
]
all_ok = True
for label, solved, expected in checks:
    rel_err = abs(solved - expected) / abs(expected)
    ok = rel_err < tol
    all_ok &= ok
    print(f"  {label:<8} rel. err {rel_err:.2e}  {'PASS' if ok else 'FAIL'}")

print(f"\n{'PASS' if all_ok else 'FAIL'}: ThermoWave's SimpleCompressor reproduces the "
      f"design-point thermodynamics implied by NASA Rotor 37's published "
      f"(PR, eta_s) to within {tol:.0e} relative error.")

assert result.converged, "network did not converge"
assert all_ok, "solved state does not match the hand-calculated design point"
