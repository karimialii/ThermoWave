"""IdealGasFluid: analytic constant-cp ideal-gas property model.

Run: .venv/bin/python examples/01_ideal_gas_fluid.py
"""

from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)

P = 101325.0  # Pa
T = 300.0  # K

h = air.enthalpy_pt(P, T)
T_back = air.temperature_ph(P, h)
rho = air.density_ph(P, h)
cp = air.cp(P, T)

print(f"fluid            : {air.name}")
print(f"P                : {P:.1f} Pa")
print(f"T                : {T:.2f} K")
print(f"h = cp*T         : {h:.2f} J/kg")
print(f"T round-tripped  : {T_back:.2f} K")
print(f"rho = P/(R*T)    : {rho:.4f} kg/m^3")
print(f"cp               : {cp:.1f} J/(kg*K)")
