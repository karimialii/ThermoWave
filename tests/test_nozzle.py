import math

import pytest

from thermowave.components.nozzle import Nozzle
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)
GAMMA = 1005.0 / (1005.0 - 287.05)
PR_CRIT = (2.0 / (GAMMA + 1.0)) ** (GAMMA / (GAMMA - 1.0))


class _FakeState:
    def __init__(self, fluid, node_values, mdots):
        self.fluid = fluid
        self._node_values = node_values
        self._mdots = mdots

    def node(self, name: str) -> tuple[float, float]:
        return self._node_values[name]

    def fluid_at(self, name: str):
        return self.fluid

    def mdot(self, name: str) -> float:
        return self._mdots[name]


def test_nozzle_rejects_eta_out_of_range():
    with pytest.raises(ValueError, match="eta"):
        Nozzle(name="n1", D=0.02, eta=0.0)
    with pytest.raises(ValueError, match="eta"):
        Nozzle(name="n1", D=0.02, eta=1.5)


def test_nozzle_rejects_gamma_not_greater_than_one():
    with pytest.raises(ValueError, match="gamma"):
        Nozzle(name="n1", D=0.02, gamma=1.0)


def test_nozzle_ports_returns_inlet_and_outlet_derived_from_name():
    noz = Nozzle(name="n1", D=0.02)
    assert noz.ports() == {"in": "n1.in", "out": "n1.out"}


def test_nozzle_subsonic_flow_matches_hand_calc():
    noz = Nozzle(name="n1", D=0.02, eta=1.0)
    P_in, T_in = 200000.0, 500.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    PR = 0.9  # well above PR_crit -> unchoked
    P_out = PR * P_in

    T_out_isentropic = T_in * PR ** ((GAMMA - 1.0) / GAMMA)
    h_out_isentropic = AIR.enthalpy_pt(P_out, T_out_isentropic)
    dh_actual = h_in - h_out_isentropic  # eta=1.0
    V = math.sqrt(2.0 * dh_actual)
    h_out_actual = h_in - dh_actual
    rho_out = AIR.density_ph(P_out, h_out_actual)
    mdot_computed = rho_out * (math.pi * 0.02**2 / 4) * V

    state = _FakeState(
        fluid=AIR,
        node_values={"n1.in": (P_in, h_in), "n1.out": (P_out, h_out_actual)},
        mdots={"n1.in": mdot_computed, "n1.out": mdot_computed},
    )
    mass_residual, energy_residual, mdot_out_residual = noz.residuals(state)
    assert math.isclose(mass_residual, 0.0, abs_tol=1e-9)
    assert math.isclose(energy_residual, 0.0, abs_tol=1e-6)
    assert math.isclose(mdot_out_residual, 0.0, abs_tol=1e-9)


def test_nozzle_report_metrics_mach_is_one_at_pr_crit():
    noz = Nozzle(name="n1", D=0.02, eta=1.0)
    P_in, T_in = 200000.0, 500.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    P_out = PR_CRIT * P_in  # exactly at the critical ratio

    state = _FakeState(
        fluid=AIR,
        node_values={"n1.in": (P_in, h_in), "n1.out": (P_out, h_in)},
        mdots={"n1.in": 1.0, "n1.out": 1.0},
    )
    metrics = noz.report_metrics(state)
    assert math.isclose(metrics["Mach [-]"], 1.0, rel_tol=1e-6)


def test_nozzle_end_to_end_subsonic_converges():
    src = Source(name="src", P=200000.0, T=500.0, mdot=None, mdot_guess=0.05)
    noz = Nozzle(name="noz", D=0.02, eta=1.0)
    snk = Sink(name="snk", P=180000.0)  # PR = 0.9, unchoked

    network = Network(fluid=AIR)
    for component in (src, noz, snk):
        network.add_component(component)
    network.connect(src, "out", noz, "in")
    network.connect(noz, "out", snk, "in")

    result = network.solve(tol=1e-9, max_iter=200, verbose=False)
    assert result.converged
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params,
    )
    metrics = noz.report_metrics(state)
    assert metrics["choked [-]"] == 0.0
    assert 0.0 < metrics["Mach [-]"] < 1.0


def test_nozzle_mass_flow_is_capped_once_choked():
    def run(P_sink):
        src = Source(name="src", P=200000.0, T=500.0, mdot=None, mdot_guess=0.1)
        noz = Nozzle(name="noz", D=0.02, eta=1.0)
        snk = Sink(name="snk", P=P_sink)
        network = Network(fluid=AIR)
        for component in (src, noz, snk):
            network.add_component(component)
        network.connect(src, "out", noz, "in")
        network.connect(noz, "out", snk, "in")
        result = network.solve(tol=1e-9, max_iter=200, verbose=False)
        state = NetworkState(
            fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
            node_mdot=result.node_mdot, params=result.params,
        )
        return noz.report_metrics(state)

    metrics_1 = run(101325.0)   # PR ~ 0.507, just below PR_crit -> choked
    metrics_2 = run(50000.0)    # much lower downstream P, still choked
    assert metrics_1["choked [-]"] == 1.0
    assert metrics_2["choked [-]"] == 1.0
    assert math.isclose(metrics_1["mdot [kg/s]"], metrics_2["mdot [kg/s]"], rel_tol=1e-6)
    assert math.isclose(metrics_1["Mach [-]"], 1.0, rel_tol=1e-6)


# ---------------------------------------------------------------------------
# Converging-diverging (de Laval) section: D_exit
# ---------------------------------------------------------------------------


def test_nozzle_rejects_d_exit_not_greater_than_throat_diameter():
    with pytest.raises(ValueError, match="D_exit"):
        Nozzle(name="n1", D=0.02, D_exit=0.02)
    with pytest.raises(ValueError, match="D_exit"):
        Nozzle(name="n1", D=0.02, D_exit=0.015)


def test_mach_from_area_ratio_matches_standard_gas_dynamics_table():
    # Textbook value (gamma=1.4): A/A* = 1.6875 <-> M = 2.0 exactly, the
    # standard reference check for this relation.
    M = Nozzle._mach_from_area_ratio(1.4, 1.6875)
    assert math.isclose(M, 2.0, rel_tol=1e-6)


def test_mach_from_area_ratio_of_one_is_sonic():
    M = Nozzle._mach_from_area_ratio(1.4, 1.0)
    assert math.isclose(M, 1.0, rel_tol=1e-6)


def test_nozzle_with_d_exit_but_unchoked_reports_no_exit_plane_metrics():
    # D_exit is inert while subsonic throughout -- same plain-converging
    # behavior, no exit-plane diagnostics reported.
    src = Source(name="src", P=200000.0, T=500.0, mdot=None, mdot_guess=0.05)
    noz = Nozzle(name="noz", D=0.02, D_exit=0.03, eta=1.0)
    snk = Sink(name="snk", P=180000.0)  # PR = 0.9, unchoked

    network = Network(fluid=AIR)
    for component in (src, noz, snk):
        network.add_component(component)
    network.connect(src, "out", noz, "in")
    network.connect(noz, "out", snk, "in")

    result = network.solve(tol=1e-9, max_iter=200, verbose=False)
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params,
    )
    metrics = noz.report_metrics(state)
    assert metrics["choked [-]"] == 0.0
    assert "Mach_exit [-]" not in metrics


def test_nozzle_with_d_exit_choked_reaches_supersonic_exit_mach():
    src = Source(name="src", P=800000.0, T=500.0, mdot=None, mdot_guess=0.1)
    noz = Nozzle(name="noz", D=0.02, D_exit=0.03, eta=1.0)
    snk = Sink(name="snk", P=50000.0)  # far below critical -> choked

    network = Network(fluid=AIR)
    for component in (src, noz, snk):
        network.add_component(component)
    network.connect(src, "out", noz, "in")
    network.connect(noz, "out", snk, "in")

    result = network.solve(tol=1e-9, max_iter=200, verbose=False)
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params,
    )
    metrics = noz.report_metrics(state)
    assert metrics["choked [-]"] == 1.0
    assert metrics["Mach [-]"] == pytest.approx(1.0, rel=1e-6)  # throat still sonic
    assert metrics["Mach_exit [-]"] > 1.0
    assert 0.0 < metrics["P_exit_ideal [Pa]"] < result.node_P["noz.in"]


def test_nozzle_exit_plane_state_matches_hand_calc_from_area_ratio():
    noz = Nozzle(name="noz", D=0.02, D_exit=0.03, eta=1.0)
    P_in, T_in = 800000.0, 500.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    P_out = 50000.0  # deep in choked regime

    area_ratio = noz._area_exit / noz._area
    M_exit_expected = Nozzle._mach_from_area_ratio(GAMMA, area_ratio)
    T_exit_expected = T_in / (1.0 + (GAMMA - 1.0) / 2.0 * M_exit_expected**2)
    P_exit_expected = P_in * (T_exit_expected / T_in) ** (GAMMA / (GAMMA - 1.0))
    h_exit_isentropic = AIR.enthalpy_pt(P_exit_expected, T_exit_expected)
    dh_exit_expected = h_in - h_exit_isentropic  # eta=1.0
    V_exit_expected = math.sqrt(2.0 * dh_exit_expected)

    state = _FakeState(
        fluid=AIR,
        node_values={"noz.in": (P_in, h_in), "noz.out": (P_out, h_in)},
        mdots={"noz.in": 1.0, "noz.out": 1.0},
    )
    flow = noz._flow(state)
    assert math.isclose(flow["P_exit_ideal"], P_exit_expected, rel_tol=1e-9)
    assert math.isclose(flow["V_exit"], V_exit_expected, rel_tol=1e-6)


def test_nozzle_mdot_is_unaffected_by_d_exit_when_choked():
    # mdot is set entirely by throat continuity -- a wider or narrower
    # diverging section must not change how much mass gets through.
    def run(D_exit):
        src = Source(name="src", P=800000.0, T=500.0, mdot=None, mdot_guess=0.1)
        noz = Nozzle(name="noz", D=0.02, D_exit=D_exit, eta=1.0)
        snk = Sink(name="snk", P=50000.0)
        network = Network(fluid=AIR)
        for component in (src, noz, snk):
            network.add_component(component)
        network.connect(src, "out", noz, "in")
        network.connect(noz, "out", snk, "in")
        result = network.solve(tol=1e-9, max_iter=200, verbose=False)
        state = NetworkState(
            fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
            node_mdot=result.node_mdot, params=result.params,
        )
        return noz.report_metrics(state)

    metrics_narrow = run(0.025)
    metrics_wide = run(0.04)
    assert math.isclose(
        metrics_narrow["mdot [kg/s]"], metrics_wide["mdot [kg/s]"], rel_tol=1e-6
    )
    # But the wider exit reaches a higher supersonic Mach (larger area ratio).
    assert metrics_wide["Mach_exit [-]"] > metrics_narrow["Mach_exit [-]"]
