import math

import pytest

from thermowave.components.pipe import Pipe
from thermowave.components.simple_compressor import SimpleCompressor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.valve import Valve
from thermowave.core.exceptions import NetworkTopologyError
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def _hand_calc_darcy_dp(fluid, L, D, f, mdot, P_in, T_in) -> float:
    h_in = fluid.enthalpy_pt(P_in, T_in)
    rho = fluid.density_ph(P_in, h_in)
    area = math.pi * D**2 / 4
    v = mdot / (rho * area)
    return f * (L / D) * (rho * v**2 / 2)


def _hand_calc_valve_dp(fluid, D, K, opening, mdot, P_in, T_in) -> float:
    h_in = fluid.enthalpy_pt(P_in, T_in)
    rho = fluid.density_ph(P_in, h_in)
    area = math.pi * D**2 / 4
    v = mdot / (rho * area)
    k_eff = K / opening**2
    return k_eff * (rho * v**2 / 2)


def _hand_calc_compressor_dh(fluid, PR, eta_s, gamma, P_in, T_in) -> float:
    h_in = fluid.enthalpy_pt(P_in, T_in)
    T_out_isentropic = T_in * PR ** ((gamma - 1.0) / gamma)
    h_out_isentropic = fluid.enthalpy_pt(P_in * PR, T_out_isentropic)
    return (h_out_isentropic - h_in) / eta_s


def test_source_pipe_sink_single_element_matches_hand_calculated_pressure_drop():
    L, D, f, mdot = 5.0, 0.2, 0.02, 1.0
    P_in, T_in = 101325.0, 300.0
    expected_dp = _hand_calc_darcy_dp(AIR, L, D, f, mdot, P_in, T_in)

    network = Network(fluid=AIR)
    src = Source(name="src", P=P_in, T=T_in, mdot=mdot)
    pipe = Pipe(name="p1", L=L, D=D, f=f, n_elem=1)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(pipe)
    network.add_component(snk)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=50)

    assert result.converged
    assert result.node_P["src.out"] == pytest.approx(P_in, rel=1e-9)
    assert result.node_P["p1.out"] == pytest.approx(P_in - expected_dp, rel=1e-4)
    assert result.node_h["p1.out"] == pytest.approx(
        result.node_h["src.out"], rel=1e-9
    )  # adiabatic


def test_source_pipe_sink_multi_element_matches_single_element_total_drop():
    L, D, f, mdot = 6.0, 0.2, 0.02, 1.0
    P_in, T_in = 101325.0, 300.0
    expected_dp = _hand_calc_darcy_dp(AIR, L, D, f, mdot, P_in, T_in)

    network = Network(fluid=AIR)
    src = Source(name="src", P=P_in, T=T_in, mdot=mdot)
    pipe = Pipe(name="p1", L=L, D=D, f=f, n_elem=3)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(pipe)
    network.add_component(snk)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=50)

    assert result.converged
    assert result.node_P["p1.out"] == pytest.approx(P_in - expected_dp, rel=1e-3)


def test_source_pipe_sink_with_heat_loss_matches_energy_balance():
    L, D, f, mdot = 5.0, 0.2, 0.02, 1.0
    Q = 2000.0  # W
    P_in, T_in = 101325.0, 300.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    expected_h_out = h_in - Q / mdot

    network = Network(fluid=AIR)
    src = Source(name="src", P=P_in, T=T_in, mdot=mdot)
    pipe = Pipe(name="p1", L=L, D=D, f=f, n_elem=1, heat_loss=Q)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(pipe)
    network.add_component(snk)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=50)

    assert result.converged
    assert result.node_h["p1.out"] == pytest.approx(expected_h_out, rel=1e-6)


def test_network_without_source_raises_topology_error_on_solve():
    network = Network(fluid=AIR)
    pipe = Pipe(name="p1", L=5.0, D=0.2, f=0.02)
    snk = Sink(name="snk")
    network.add_component(pipe)
    network.add_component(snk)
    network.connect(pipe, "out", snk, "in")

    with pytest.raises(NetworkTopologyError):
        network.solve()


def test_source_valve_sink_matches_hand_calculated_pressure_drop():
    D, K, mdot = 0.1, 5.0, 0.5
    P_in, T_in = 101325.0, 300.0
    expected_dp = _hand_calc_valve_dp(AIR, D, K, 1.0, mdot, P_in, T_in)

    network = Network(fluid=AIR)
    src = Source(name="src", P=P_in, T=T_in, mdot=mdot)
    valve = Valve(name="v1", D=D, K=K)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(valve)
    network.add_component(snk)
    network.connect(src, "out", valve, "in")
    network.connect(valve, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=50)

    assert result.converged
    assert result.node_P["v1.out"] == pytest.approx(P_in - expected_dp, rel=1e-4)
    assert result.node_h["v1.out"] == pytest.approx(
        result.node_h["src.out"], rel=1e-9
    )  # isenthalpic


def test_source_pipe_valve_sink_partially_closed_matches_hand_calc():
    L, D_pipe, f, mdot = 5.0, 0.2, 0.02, 0.3
    D_valve, K, opening = 0.15, 2.0, 0.6
    P_in, T_in = 101325.0, 300.0
    expected_pipe_dp = _hand_calc_darcy_dp(AIR, L, D_pipe, f, mdot, P_in, T_in)

    network = Network(fluid=AIR)
    src = Source(name="src", P=P_in, T=T_in, mdot=mdot)
    pipe = Pipe(name="p1", L=L, D=D_pipe, f=f)
    valve = Valve(name="v1", D=D_valve, K=K, opening=opening)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(pipe)
    network.add_component(valve)
    network.add_component(snk)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", valve, "in")
    network.connect(valve, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=50)
    assert result.converged

    P_after_pipe = P_in - expected_pipe_dp
    expected_valve_dp = _hand_calc_valve_dp(AIR, D_valve, K, opening, mdot, P_after_pipe, T_in)
    assert result.node_P["p1.out"] == pytest.approx(P_after_pipe, rel=1e-3)
    assert result.node_P["v1.out"] == pytest.approx(
        P_after_pipe - expected_valve_dp, rel=1e-3
    )


def test_source_simple_compressor_sink_matches_isentropic_efficiency_hand_calc():
    PR, eta_s, gamma = 3.0, 0.8, 1005.0 / (1005.0 - 287.05)
    mdot = 1.0
    P_in, T_in = 101325.0, 300.0
    expected_dh = _hand_calc_compressor_dh(AIR, PR, eta_s, gamma, P_in, T_in)
    h_in = AIR.enthalpy_pt(P_in, T_in)

    network = Network(fluid=AIR)
    src = Source(name="src", P=P_in, T=T_in, mdot=mdot)
    comp = SimpleCompressor(name="c1", PR=PR, eta_s=eta_s, gamma=gamma)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(comp)
    network.add_component(snk)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=50)

    assert result.converged
    assert result.node_P["c1.out"] == pytest.approx(PR * P_in, rel=1e-9)
    assert result.node_h["c1.out"] == pytest.approx(h_in + expected_dh, rel=1e-4)
