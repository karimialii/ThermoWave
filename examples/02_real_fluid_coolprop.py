"""CoolPropFluid: real-fluid properties backed by CoolProp's PropsSI.

Requires the optional 'coolprop' extra: pip install thermowave[coolprop]
Run: .venv/bin/python examples/02_real_fluid_coolprop.py
"""

try:
    from thermowave.fluids.real_fluid import CoolPropFluid
except ImportError as exc:
    raise SystemExit(f"Skipping example, CoolProp not installed: {exc}")

water = CoolPropFluid(name="Water")

P = 1.0e5  # Pa, 1 bar
T = 473.15  # K, superheated steam at 1 bar

h = water.enthalpy_pt(P, T)
T_back = water.temperature_ph(P, h)
rho = water.density_ph(P, h)

print(f"fluid                    : {water.name}")
print(f"P                        : {P:.1f} Pa")
print(f"T                        : {T:.2f} K")
print(f"h (PropsSI)              : {h:.2f} J/kg")
print(f"T round-tripped from h   : {T_back:.2f} K")
print(f"rho (PropsSI, from P,h)  : {rho:.4f} kg/m^3")

clamped_low = water._clamp_pressure(1.0)
clamped_high = water._clamp_pressure(1.0e9)
print(f"pressure clamp floor     : {clamped_low:.1f} Pa (P_min={water.P_min:.1e})")
print(f"pressure clamp ceiling   : {clamped_high:.1f} Pa (P_max={water.P_max:.1e})")
