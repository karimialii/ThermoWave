"""Source: boundary component fixing outlet (P, T) and the network's mass flow rate.

A component's ports() are auto-derived from its name (e.g. "src1.out");
Network.connect() is what wires them to other components (see example 07).

Run: .venv/bin/python examples/04_source_component.py
"""

from thermowave.components.source import Source
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)

src = Source(name="src1", P=101325.0, T=300.0, mdot=1.0)

print(f"ports              : {src.ports()}")
print(f"fixed_node_mdot()  : {src.fixed_node_mdot()} kg/s")

fixed = src.fixed_node_values(air)
node_name, (P, h) = next(iter(fixed.items()))
print(f"fixed_node_values: node {node_name!r} -> P={P:.1f} Pa, h={h:.1f} J/kg")

# Source also accepts non-SI input via the settings singleton (see example 03);
# construction-time conversion means fixed_node_values() always returns SI.
