import math

import pytest

pytest.importorskip("CoolProp")

from thermowave.components.sink import Sink  # noqa: E402
from thermowave.components.source import Source  # noqa: E402
from thermowave.components.steam_turbine import SteamTurbine  # noqa: E402
from thermowave.core.network import Network, NetworkState  # noqa: E402
from thermowave.fluids.real_fluid import CoolPropFluid  # noqa: E402

WATER = CoolPropFluid(name="Water")


class _FakeState:
    def __init__(self, fluid, node_values, mdots):
        self.fluid = fluid
        self._node_values = node_values
        self._mdots = mdots

    def node(self, name):
        return self._node_values[name]

    def fluid_at(self, name):
        return self.fluid

    def mdot(self, name):
        return self._mdots[name]


def test_rejects_both_or_neither_of_p_out_and_pr():
    with pytest.raises(ValueError, match="exactly one"):
        SteamTurbine(name="st", P_out=1.0e4, PR=100.0)
    with pytest.raises(ValueError, match="exactly one"):
        SteamTurbine(name="st")


def test_isentropic_drop_hand_calc_with_efficiency():
    eta_s = 0.85
    st = SteamTurbine(name="st", P_out=1.0e4, eta_s=eta_s)
    P_in = 1.0e6
    h_in = WATER.enthalpy_pt(P_in, 600.0)  # superheated
    P_out = 1.0e4
    s_in = WATER.entropy_ph(P_in, h_in)
    h_out_isentropic = WATER.enthalpy_ps(P_out, s_in)
    h_out = h_in - eta_s * (h_in - h_out_isentropic)
    state = _FakeState(
        fluid=WATER,
        node_values={"st.in": (P_in, h_in), "st.out": (P_out, h_out)},
        mdots={"st.in": 1.0, "st.out": 1.0},
    )
    momentum, energy, mass = st.residuals(state)
    assert math.isclose(momentum, 0.0, abs_tol=1e-3)
    assert math.isclose(energy, 0.0, abs_tol=1e-6)
    assert math.isclose(mass, 0.0, abs_tol=1e-12)


def test_lower_efficiency_extracts_less_work():
    P_in, P_out = 1.0e6, 1.0e4
    h_in = WATER.enthalpy_pt(P_in, 600.0)
    s_in = WATER.entropy_ph(P_in, h_in)
    h_out_isentropic = WATER.enthalpy_ps(P_out, s_in)

    def work(eta_s):
        st = SteamTurbine(name="st", P_out=P_out, eta_s=eta_s)
        h_out = h_in - eta_s * (h_in - h_out_isentropic)
        state = _FakeState(
            fluid=WATER,
            node_values={"st.in": (P_in, h_in), "st.out": (P_out, h_out)},
            mdots={"st.in": 1.0, "st.out": 1.0},
        )
        return st.report_metrics(state)["power [W]"]

    assert work(0.7) < work(0.9)


def test_end_to_end_expands_into_wet_region_and_reports_quality():
    src = Source(name="src", P=1.0e6, T=600.0, mdot=3.0)  # superheated steam
    st = SteamTurbine(name="st", P_out=1.0e4, eta_s=0.85)
    snk = Sink(name="snk")

    network = Network(fluid=WATER)
    for c in (src, st, snk):
        network.add_component(c)
    network.connect(src, "out", st, "in")
    network.connect(st, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=100)
    assert result.converged
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params,
    )
    metrics = st.report_metrics(state)
    assert metrics["power [W]"] > 0.0
    # Expansion from superheated 10 bar to 0.1 bar lands in the wet region.
    assert 0.0 < metrics["x_out [-]"] < 1.0
