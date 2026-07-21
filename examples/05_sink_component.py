"""Sink: open boundary terminating a network branch. Fixes no state itself.

Run: .venv/bin/python examples/05_sink_component.py
"""

from thermowave.components.sink import Sink
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)

snk = Sink(name="snk1")

print(f"ports              : {snk.ports()}")
print(f"fixed_node_values()  : {snk.fixed_node_values(air)}  (empty — Sink fixes nothing)")
print(f"fixed_node_mdot()    : {snk.fixed_node_mdot()}  (empty — Sink doesn't fix flow rate)")
