"""Valve: K-factor pressure drop, isenthalpic (no heat/work), with an optional
opening fraction (0, 1] that scales resistance up as the valve closes.

Run: .venv/bin/python examples/07_valve_component.py
"""

import math

from thermowave.components.valve import Valve
from thermowave.core.network import NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)

D, K, mdot = 0.1, 5.0, 0.5
P_in, T_in = 101325.0, 300.0
h_in = air.enthalpy_pt(P_in, T_in)

valve = Valve(name="v1", D=D, K=K)
print(f"ports: {valve.ports()}")

inlet_node = valve.ports()["in"]
outlet_node = valve.ports()["out"]

area = math.pi * D**2 / 4
rho = air.density_ph(P_in, h_in)
v = mdot / (rho * area)
dp_fully_open = K * (rho * v**2 / 2)

state = NetworkState(
    fluid=air,
    node_P={inlet_node: P_in, outlet_node: P_in - dp_fully_open},
    node_h={inlet_node: h_in, outlet_node: h_in},  # isenthalpic
    node_mdot={inlet_node: mdot, outlet_node: mdot},
)
momentum_residual, energy_residual, _mass_residual = valve.residuals(state)
print(f"fully open (opening=1.0) dp guess: {dp_fully_open:.2f} Pa")
print(f"momentum residual: {momentum_residual:.3e}  (~0 means the guess is correct)")
print(f"energy residual  : {energy_residual:.3e}")

# Closing the valve raises the effective loss coefficient (K / opening**2),
# so the same guessed dp is no longer a solution — residual moves off zero.
half_closed = Valve(name="v1", D=D, K=K, opening=0.5)
_momentum_residual_half, _, _ = half_closed.residuals(state)
print(f"\nsame guess through a half-closed valve (opening=0.5):")
print(f"momentum residual: {_momentum_residual_half:.3e}  (nonzero: needs more dp now)")
