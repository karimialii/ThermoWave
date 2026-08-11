import math

from thermowave.components.pipe import Pipe
from thermowave.components.simple_heater import SimpleHeater
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.constants import STANDARD_GRAVITY_M_S2
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def _build_loop(H=5.0, T_cold=300.0, T_hot=340.0):
    """Source(mdot=None)/Sink pair at the same P (a pressure-fixing
    reservoir, the way a real closed thermosiphon needs a pressurizer or
    expansion tank to anchor absolute pressure) -> heater -> heated riser
    (z: 0 -> H) -> cooler -> cooled downcomer (z: H -> 0) -> Sink.

    No pump anywhere: mdot is left free (Source(mdot=None)) and must come
    entirely out of the loop's own momentum balance -- friction (which
    needs mdot to be nonzero to resist) against the buoyancy created by the
    riser running hot/light and the downcomer running cold/dense.
    """
    net = Network(fluid=AIR)
    src = Source(name="src", P=101_325.0, T=T_cold, mdot=None, mdot_guess=0.05)
    heater = SimpleHeater(name="heater", T_out=T_hot)
    riser = Pipe(name="riser", L=H, D=0.15, roughness=1e-5, mu=1.8e-5, z_in=0.0, z_out=H)
    cooler = SimpleHeater(name="cooler", T_out=T_cold)
    downcomer = Pipe(name="downcomer", L=H, D=0.15, roughness=1e-5, mu=1.8e-5, z_in=H, z_out=0.0)
    snk = Sink(name="snk", P=101_325.0)

    for c in (src, heater, riser, cooler, downcomer, snk):
        net.add_component(c)
    net.connect(src, "out", heater, "in")
    net.connect(heater, "out", riser, "in")
    net.connect(riser, "out", cooler, "in")
    net.connect(cooler, "out", downcomer, "in")
    net.connect(downcomer, "out", snk, "in")
    return net, dict(src=src, heater=heater, riser=riser, cooler=cooler, downcomer=downcomer, snk=snk)


def test_natural_circulation_loop_converges_to_nonzero_mdot():
    net, c = _build_loop()
    assert net.check_wiring() == []

    result = net.solve(tol=1e-9, max_iter=200, verbose=False, progress=False)
    assert result.converged

    mdot = result.node_mdot["riser.in"]
    assert mdot > 1.0e-4  # genuinely circulating, not stuck at a trivial root


def test_larger_temperature_difference_drives_more_flow():
    net_small, _ = _build_loop(T_cold=300.0, T_hot=310.0)
    net_large, _ = _build_loop(T_cold=300.0, T_hot=360.0)

    result_small = net_small.solve(tol=1e-9, max_iter=200, verbose=False, progress=False)
    result_large = net_large.solve(tol=1e-9, max_iter=200, verbose=False, progress=False)
    assert result_small.converged and result_large.converged

    mdot_small = result_small.node_mdot["riser.in"]
    mdot_large = result_large.node_mdot["riser.in"]
    assert mdot_large > mdot_small > 0.0


def test_circulation_matches_analytic_buoyancy_vs_friction_balance():
    H = 5.0
    net, c = _build_loop(H=H, T_cold=300.0, T_hot=340.0)
    result = net.solve(tol=1e-9, max_iter=300, verbose=False, progress=False)
    assert result.converged

    state = result.state()
    mdot = result.node_mdot["riser.in"]

    rho_hot = AIR.density_ph(*state.node("riser.in"))
    rho_cold = AIR.density_ph(*state.node("downcomer.in"))
    area = math.pi * 0.15**2 / 4

    def dp_friction(pipe, rho):
        v = mdot / (rho * area)
        f = pipe._friction_factor(rho, v)
        return f * (pipe.L / pipe.D) * (rho * v**2 / 2)

    lhs = dp_friction(c["riser"], rho_hot) + dp_friction(c["downcomer"], rho_cold)
    rhs = STANDARD_GRAVITY_M_S2 * H * (rho_cold - rho_hot)
    assert math.isclose(lhs, rhs, rel_tol=1e-4)
