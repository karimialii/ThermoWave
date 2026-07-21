"""Compressor: every supported known/unknown combination, driven by the same
T100 Comp.cop map, to show the solver finding N regardless of which quantity
you actually know.

  1. N given directly (as always).
  2. N unknown, target PR       -> Setpoint on "PR [-]"
  3. N unknown, target power    -> Setpoint on "power [W]"
  4. N unknown, target eta_s    -> Setpoint on "eta_s [-]"

Case 4 is numerically harder and included deliberately: compressor
efficiency vs speed is a peaked (non-monotonic) curve on this map, so near
the peak the Jacobian is close to singular and the default damping=1.0
overshoots. Lower damping (slower, more steps, but robust) is needed — a
concrete example of why "give me any known, I'll solve the rest" isn't free:
the equation you get has to actually be well-posed for what you removed.

Run: .venv/bin/python examples/14_compressor_all_cases.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.setpoint import Setpoint
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)


def run(label, N, target_metric, value, damping=1.0, max_iter=100):
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    comp = Compressor(name="c1", map_path="T100 Comp.cop", gamma=gamma, N=N)
    snk = Sink(name="snk")

    network = Network(fluid=air)
    network.add_component(src)
    network.add_component(comp)
    if target_metric is not None:
        sp = Setpoint(
            name="sp1", component=comp, free_param="N", target_metric=target_metric, value=value
        )
        network.add_component(sp)
    network.add_component(snk)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", snk, "in")

    print(f"\n--- {label} ---")
    result = network.solve(tol=1e-8, max_iter=max_iter, damping=damping, verbose=True)
    result.print_report()


run("1. N = 65000 rev/min (given directly)", N=65000.0, target_metric=None, value=None)
run("2. target PR = 3.8", N=None, target_metric="PR [-]", value=3.909)
run("3. target power = 112 kW", N=None, target_metric="power [W]", value=65000.0)
run(
    "4. target eta_s = 0.70 (harder: non-monotonic curve, needs damping)",
    N=None,
    target_metric="eta_s [-]",
    value=0.70,
    damping=0.3,
    max_iter=200,
)
