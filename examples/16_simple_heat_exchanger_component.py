"""SimpleHeatExchanger: 0D two-stream heat exchanger, fixed effectiveness.

Hot air stream (Source -> hot_in/hot_out -> Sink) exchanges heat with a cold
air stream (Source -> cold_in/cold_out -> Sink) at a fixed effectiveness
(given directly, not derived from UA/NTU/flow arrangement), with each stream
also seeing its own K-factor pressure drop. Both streams happen to use the
same fluid model here (air), but they're wired as two independent flow paths
through the same component, each with its own mdot.

Run: .venv/bin/python examples/16_simple_heat_exchanger_component.py
"""

from thermowave.components.simple_heat_exchanger import SimpleHeatExchanger
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)

hot_src = Source(name="hot_src", P=300000.0, T=500.0, mdot=1.0)
cold_src = Source(name="cold_src", P=300000.0, T=300.0, mdot=1.2)
hx = SimpleHeatExchanger(
    name="hx", effectiveness=0.7, PR_cold=0.95, PR_hot=1.0
)
hot_snk = Sink(name="hot_snk")
cold_snk = Sink(name="cold_snk")

network = Network(fluid=air)
for component in (hot_src, cold_src, hx, hot_snk, cold_snk):
    network.add_component(component)

network.connect(hot_src, "out", hx, "hot_in")
network.connect(hx, "hot_out", hot_snk, "in")
network.connect(cold_src, "out", hx, "cold_in")
network.connect(hx, "cold_out", cold_snk, "in")

result = network.solve(tol=1e-8, max_iter=50, verbose=True)
result.print_report()

state = NetworkState(
    fluid=air,
    node_P=result.node_P,
    node_h=result.node_h,
    node_mdot=result.node_mdot,
    params=result.params,
)
metrics = hx.report_metrics(state)
print(f"\nHeat duty Q = {metrics['power [W]']:.1f} W")
print(f"T_hot_in = {metrics['T_hot_in [K]']:.2f} K, T_cold_in = {metrics['T_cold_in [K]']:.2f} K")
