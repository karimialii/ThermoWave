"""Pipe: Darcy-Weisbach friction pressure drop + optional heat loss, discretized
into n_elem sub-elements. Each element contributes a momentum and an energy
residual, evaluated against a NetworkState of node (P, h) values.

Run: .venv/bin/python examples/06_pipe_component.py
"""

import math

from thermowave.components.pipe import Pipe
from thermowave.core.network import NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)

L, D, f, mdot = 5.0, 0.2, 0.02, 1.0
P_in, T_in = 101325.0, 300.0
h_in = air.enthalpy_pt(P_in, T_in)

pipe = Pipe(name="p1", L=L, D=D, f=f, n_elem=1)
print(f"ports            : {pipe.ports()}")
print(f"internal_nodes() : {pipe.internal_nodes()}  (empty for a single element)")

# Guess an outlet state and check the residuals (should be ~0 at the true answer,
# which is exactly what the Newton solver drives toward inside Network.solve()).
inlet_node = pipe.ports()["in"]
outlet_node = pipe.ports()["out"]

area = math.pi * D**2 / 4
rho = air.density_ph(P_in, h_in)
v = mdot / (rho * area)
dp_guess = f * (L / D) * (rho * v**2 / 2)

state = NetworkState(
    fluid=air,
    node_P={inlet_node: P_in, outlet_node: P_in - dp_guess},
    node_h={inlet_node: h_in, outlet_node: h_in},  # adiabatic: no heat loss
    node_mdot={inlet_node: mdot, outlet_node: mdot},
)
momentum_residual, energy_residual, _mass_residual = pipe.residuals(state)
print(f"friction dp guess: {dp_guess:.4f} Pa")
print(f"momentum residual: {momentum_residual:.3e}  (~0 means the guess is correct)")
print(f"energy residual  : {energy_residual:.3e}")

# A 3-element pipe exposes internal nodes for the mid-element states.
multi_pipe = Pipe(name="p2", L=6.0, D=D, f=f, n_elem=3)
print(f"\n3-element internal_nodes(): {multi_pipe.internal_nodes()}")
