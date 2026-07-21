"""SimpleCompressor: fixed pressure ratio + isentropic efficiency, no map.

Uses the ideal-gas isentropic relation T2s = T1 * PR**((gamma-1)/gamma), so
gamma is a required input (not derived from the fluid model) — this keeps it
usable with any BaseFluid, since real-fluid models don't expose entropy or
gamma directly. A map-based Compressor (Flownex-style iso-speed pressure
ratio / efficiency vs non-dimensional mass flow curves) is planned separately.

Run: .venv/bin/python examples/09_simple_compressor_component.py
"""

from thermowave.components.simple_compressor import SimpleCompressor
from thermowave.core.network import NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)  # cp / (cp - R), matches real air (~1.4)

PR, eta_s = 3.0, 0.8
P_in, T_in = 101325.0, 300.0
h_in = air.enthalpy_pt(P_in, T_in)

comp = SimpleCompressor(name="c1", PR=PR, eta_s=eta_s, gamma=gamma)
print(f"ports: {comp.ports()}")

inlet_node = comp.ports()["in"]
outlet_node = comp.ports()["out"]

T_out_isentropic = T_in * PR ** ((gamma - 1.0) / gamma)
h_out_isentropic = air.enthalpy_pt(P_in * PR, T_out_isentropic)
dh_actual = (h_out_isentropic - h_in) / eta_s
h_out_guess = h_in + dh_actual
P_out_guess = PR * P_in

state = NetworkState(
    fluid=air,
    node_P={inlet_node: P_in, outlet_node: P_out_guess},
    node_h={inlet_node: h_in, outlet_node: h_out_guess},
    node_mdot={inlet_node: 1.0, outlet_node: 1.0},
)
momentum_residual, energy_residual, _mass_residual = comp.residuals(state)

print(f"P_in                    : {P_in:.1f} Pa")
print(f"P_out = PR * P_in       : {P_out_guess:.1f} Pa")
print(f"T_out isentropic        : {T_out_isentropic:.2f} K")
print(f"dh_actual = dh_s/eta_s  : {dh_actual:.2f} J/kg")
print(f"momentum residual       : {momentum_residual:.3e}  (~0 means the guess is correct)")
print(f"energy residual         : {energy_residual:.3e}")

power = 1.0 * dh_actual  # mdot * dh_actual
print(f"\nshaft power at mdot=1.0 kg/s: {power / 1000:.2f} kW")
