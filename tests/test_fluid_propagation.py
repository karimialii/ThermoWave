import math

import pytest

from thermowave.components.base_component import BaseComponent
from thermowave.components.pipe import Pipe
from thermowave.components.sensor import Sensor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.base_fluid import BaseFluid
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)
OTHER = IdealGasFluid(name="other", R=300.0, cp=1100.0)


# ---------------------------------------------------------------------------
# NetworkState.fluid_at()
# ---------------------------------------------------------------------------


def test_fluid_at_falls_back_to_default_fluid_when_node_has_no_override():
    state = NetworkState(fluid=AIR, node_P={}, node_h={}, node_mdot={})
    assert state.fluid_at("anything.in") is AIR


def test_fluid_at_returns_node_specific_override_when_present():
    state = NetworkState(
        fluid=AIR, node_P={}, node_h={}, node_mdot={}, node_fluid={"x.out": OTHER}
    )
    assert state.fluid_at("x.out") is OTHER
    assert state.fluid_at("x.in") is AIR


def test_base_component_outlet_fluid_default_is_pass_through():
    class _Dummy(BaseComponent):
        name = "d1"

        def ports(self):
            return {"in": "d1.in", "out": "d1.out"}

        def residuals(self, state):
            return []

    state = NetworkState(fluid=AIR, node_P={}, node_h={}, node_mdot={})
    assert _Dummy().outlet_fluid(state, ("in", "out"), AIR) is None


# ---------------------------------------------------------------------------
# Network._resolve_node_fluid(): propagation mechanics
# ---------------------------------------------------------------------------


class _CompositionSwap(BaseComponent):
    """Test-only component that swaps in a different BaseFluid on its
    outlet, unconditionally -- stands in for Combustor without needing
    cantera, to test the propagation loop itself in isolation."""

    def __init__(self, name, new_fluid):
        self.name = name
        self.new_fluid = new_fluid
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"

    def ports(self):
        return {"in": self._inlet_node, "out": self._outlet_node}

    def residuals(self, state):
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        mass_residual = state.mdot(self._outlet_node) - state.mdot(self._inlet_node)
        return [P_out - P_in, h_out - h_in, mass_residual]

    def outlet_fluid(self, state, pair, inlet_fluid):
        if pair != ("in", "out"):
            return None
        return self.new_fluid


def _chain_network(register_order):
    src = Source(name="src", P=200000.0, T=400.0, mdot=1.0)
    swap = _CompositionSwap(name="swap", new_fluid=OTHER)
    pipe = Pipe(name="pipe", L=1.0, D=0.1, f=0.02)
    snk = Sink(name="snk")

    network = Network(fluid=AIR)
    components = {"src": src, "swap": swap, "pipe": pipe, "snk": snk}
    for key in register_order:
        network.add_component(components[key])
    network.connect(src, "out", swap, "in")
    network.connect(swap, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")
    return network, src, swap, pipe, snk


def test_resolve_node_fluid_propagates_pass_through_when_nothing_changes_it():
    src = Source(name="src", P=200000.0, T=400.0, mdot=1.0)
    pipe = Pipe(name="pipe", L=1.0, D=0.1, f=0.02)
    snk = Sink(name="snk")
    network = Network(fluid=AIR)
    for c in (src, pipe, snk):
        network.add_component(c)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")

    result = network.solve(tol=1e-9, max_iter=100)
    # Nothing in this network overrides outlet_fluid(), so every node either
    # has no node_fluid entry, or (if present) points at the same fluid.
    for node, fluid in result.node_fluid.items():
        assert fluid is AIR


def test_resolve_node_fluid_propagates_composition_change_downstream_only():
    network, src, swap, pipe, snk = _chain_network(["src", "swap", "pipe", "snk"])
    result = network.solve(tol=1e-9, max_iter=100)

    assert result.node_fluid["src.out"] is AIR
    assert result.node_fluid["swap.out"] is OTHER
    assert result.node_fluid["pipe.out"] is OTHER


def test_resolve_node_fluid_is_order_independent():
    # Components registered out of upstream-to-downstream order must still
    # resolve correctly -- the fixed-point loop (same technique as the
    # solver's own P/h warm-start propagation) is designed for exactly this.
    network, src, swap, pipe, snk = _chain_network(["snk", "pipe", "swap", "src"])
    result = network.solve(tol=1e-9, max_iter=100)

    assert result.node_fluid["swap.out"] is OTHER
    assert result.node_fluid["pipe.out"] is OTHER


def test_resolve_node_fluid_downstream_component_actually_reads_the_swap():
    # Prove this isn't just bookkeeping: Pipe's own residual (via
    # state.fluid_at()) must use OTHER's density for the segment after the
    # swap, not AIR's -- i.e. the propagated fluid actually feeds physics,
    # not just report_metrics().
    network, src, swap, pipe, snk = _chain_network(["src", "swap", "pipe", "snk"])
    result = network.solve(tol=1e-9, max_iter=100)

    P_in, h_in = result.node_P["pipe.in"], result.node_h["pipe.in"]
    rho_other = OTHER.density_ph(P_in, h_in)
    rho_air = AIR.density_ph(P_in, h_in)
    assert not math.isclose(rho_other, rho_air, rel_tol=1e-3)

    mdot = result.node_mdot["pipe.in"]
    area = math.pi * pipe.D**2 / 4
    v = mdot / (rho_other * area)
    expected_dp = pipe.f * (pipe.L / pipe.D) * (rho_other * v**2 / 2)
    assert math.isclose(
        result.node_P["pipe.in"] - result.node_P["pipe.out"], expected_dp, rel_tol=1e-6
    )


# ---------------------------------------------------------------------------
# Combustor.outlet_fluid(): gating + live composition propagation
# ---------------------------------------------------------------------------

cantera = pytest.importorskip("cantera")


def test_combustor_outlet_fluid_is_none_for_non_cantera_inlet():
    from thermowave.components.combustor import Combustor

    comb = Combustor(name="cc1", mdot_fuel=0.02)
    state = NetworkState(
        fluid=AIR,
        node_P={"cc1.in": 300000.0},
        node_h={"cc1.in": AIR.enthalpy_pt(300000.0, 500.0)},
        node_mdot={"cc1.in": 1.0},
    )
    assert comb.outlet_fluid(state, ("in", "out"), AIR) is None


def test_combustor_outlet_fluid_returns_working_fluid_for_cantera_inlet():
    from thermowave.components.combustor import Combustor
    from thermowave.fluids.cantera_fluid import CanteraFluid

    cantera_air = CanteraFluid(name="cantera_air", composition="O2:0.21, N2:0.79")
    comb = Combustor(name="cc1", mdot_fuel=0.02, fuel="CH4")
    P_in, T_in = 300000.0, 500.0
    h_in = cantera_air.enthalpy_pt(P_in, T_in)
    state = NetworkState(
        fluid=cantera_air,
        node_P={"cc1.in": P_in},
        node_h={"cc1.in": h_in},
        node_mdot={"cc1.in": 1.0},
    )

    product_fluid = comb.outlet_fluid(state, ("in", "out"), cantera_air)
    assert product_fluid is not None
    assert product_fluid is not cantera_air

    # Round-trip consistency: the product fluid's own temperature_ph/
    # enthalpy_pt must agree with each other at some arbitrary (P, T).
    P_out, T_probe = 0.96 * P_in, 1200.0
    h_probe = product_fluid.enthalpy_pt(P_out, T_probe)
    T_recovered = product_fluid.temperature_ph(P_out, h_probe)
    assert math.isclose(T_recovered, T_probe, rel_tol=1e-6)


def test_combustor_mixture_cache_avoids_redundant_equilibrate_calls():
    from thermowave.components.combustor import Combustor
    from thermowave.fluids.cantera_fluid import CanteraFluid

    cantera_air = CanteraFluid(name="cantera_air", composition="O2:0.21, N2:0.79")
    comb = Combustor(name="cc1", mdot_fuel=0.02, fuel="CH4")
    P_in, T_in = 300000.0, 500.0
    state = NetworkState(
        fluid=cantera_air,
        node_P={"cc1.in": P_in},
        node_h={"cc1.in": cantera_air.enthalpy_pt(P_in, T_in)},
        node_mdot={"cc1.in": 1.0},
    )

    calls = []
    real_equilibrate = comb._equilibrate

    def _counting_equilibrate(*args, **kwargs):
        calls.append(1)
        return real_equilibrate(*args, **kwargs)

    comb._equilibrate = _counting_equilibrate
    comb.outlet_fluid(state, ("in", "out"), cantera_air)
    comb.outlet_fluid(state, ("in", "out"), cantera_air)
    comb.product_composition(state)
    assert len(calls) == 1


def test_combustor_end_to_end_composition_propagates_to_downstream_pipe():
    from thermowave.components.combustor import Combustor
    from thermowave.fluids.cantera_fluid import CanteraFluid

    cantera_air = CanteraFluid(name="cantera_air", composition="O2:0.21, N2:0.79")
    src = Source(name="src", P=300000.0, T=500.0, mdot=1.0)
    comb = Combustor(name="cc1", PR=0.96, mdot_fuel=0.02, fuel="CH4")
    pipe = Pipe(name="pipe", L=2.0, D=0.15, f=0.02)
    snk = Sink(name="snk")

    network = Network(fluid=cantera_air)
    for c in (src, comb, pipe, snk):
        network.add_component(c)
    network.connect(src, "out", comb, "in")
    network.connect(comb, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")

    result = network.solve(tol=1e-6, max_iter=300)

    product_fluid = result.node_fluid["pipe.in"]
    assert product_fluid is not cantera_air

    # Live wiring proof: Pipe's actual pressure drop must match what the
    # *product* fluid's density predicts at its own inlet state, not what
    # plain pre-combustion air would predict at the same (P, h).
    P_in, h_in = result.node_P["pipe.in"], result.node_h["pipe.in"]
    mdot = result.node_mdot["pipe.in"]
    area = math.pi * pipe.D**2 / 4

    rho_product = product_fluid.density_ph(P_in, h_in)
    rho_air = cantera_air.density_ph(P_in, h_in)
    assert not math.isclose(rho_product, rho_air, rel_tol=1e-2)

    v_product = mdot / (rho_product * area)
    expected_dp = pipe.f * (pipe.L / pipe.D) * (rho_product * v_product**2 / 2)
    actual_dp = result.node_P["pipe.in"] - result.node_P["pipe.out"]
    assert math.isclose(actual_dp, expected_dp, rel_tol=1e-6)


def test_combustor_downstream_sensor_reads_product_fluid_temperature():
    from thermowave.components.combustor import Combustor
    from thermowave.fluids.cantera_fluid import CanteraFluid

    cantera_air = CanteraFluid(name="cantera_air", composition="O2:0.21, N2:0.79")
    src = Source(name="src", P=300000.0, T=500.0, mdot=1.0)
    comb = Combustor(name="cc1", PR=0.96, mdot_fuel=0.02, fuel="CH4")
    sensor = Sensor(name="s1")
    snk = Sink(name="snk")

    network = Network(fluid=cantera_air)
    for c in (src, comb, sensor, snk):
        network.add_component(c)
    network.connect(src, "out", comb, "in")
    network.connect(comb, "out", sensor, "tap")
    network.connect(sensor, "tap", snk, "in")

    result = network.solve(tol=1e-6, max_iter=300)
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params, node_fluid=result.node_fluid,
    )
    metrics = sensor.report_metrics(state)

    product_fluid = result.node_fluid["cc1.out"]
    expected_T = product_fluid.temperature_ph(metrics["P [Pa]"], metrics["h [J/kg]"])
    assert math.isclose(metrics["T [K]"], expected_T, rel_tol=1e-9)
    # And it must differ from reading the same (P, h) via plain air --
    # otherwise the sensor would just be silently ignoring the propagated
    # composition.
    T_via_air = cantera_air.temperature_ph(metrics["P [Pa]"], metrics["h [J/kg]"])
    assert not math.isclose(metrics["T [K]"], T_via_air, rel_tol=1e-3)
