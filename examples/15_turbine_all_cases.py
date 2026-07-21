"""Turbine: every supported known/unknown combination, driven by the same
T100 Turb.tur map, mirroring 14_compressor_all_cases.py.

  1. N given directly (as always).
  2. N unknown, target PR       -> Setpoint on "PR [-]"
  3. N unknown, target power    -> Setpoint on "power [W]"
  4. N unknown, target eta_s    -> Setpoint on "eta_s [-]"

All three target cases need reduced damping here: the initial N guess
(seeded from the map's own mid-speed) starts far enough from the converged
operating point that full Newton steps overshoot into a region of the map
where the corrected mass flow clamps flat, singularizing the Jacobian.
Lower damping trades iteration count for robustness.

Run: .venv/bin/python examples/15_turbine_all_cases.py
"""

from thermowave.components.setpoint import Setpoint
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)


def run(label, N, target_metric, value, damping=1.0, max_iter=100):
    src = Source(name="src", P=400000.0, T=1150.0, mdot=0.65)
    turb = Turbine(name="t1", map_path="T100 Turb.tur", gamma=gamma, N=N)
    snk = Sink(name="snk")

    network = Network(fluid=air)
    network.add_component(src)
    network.add_component(turb)
    if target_metric is not None:
        sp = Setpoint(
            name="sp1", component=turb, free_param="N", target_metric=target_metric, value=value
        )
        network.add_component(sp)
    network.add_component(snk)
    network.connect(src, "out", turb, "in")
    network.connect(turb, "out", snk, "in")

    print(f"\n--- {label} ---")
    result = network.solve(tol=1e-8, max_iter=max_iter, damping=damping, verbose=True)
    result.print_report()


run("1. N = 47000 rev/min (given directly)", N=47000.0, target_metric=None, value=None)
run(
    "2. target PR = 1.7 (needs damping)",
    N=None,
    target_metric="PR [-]",
    value=1.7,
    damping=0.5,
    max_iter=200,
)
run(
    "3. target power = 90 kW (needs damping)",
    N=None,
    target_metric="power [W]",
    value=90000.0,
    damping=0.5,
    max_iter=200,
)
run(
    "4. target eta_s = 0.85 (harder: needs more damping)",
    N=None,
    target_metric="eta_s [-]",
    value=0.85,
    damping=0.3,
    max_iter=300,
)
