import pytest

pytest.importorskip("CoolProp")

from thermowave.components.simple_heater import SimpleHeater  # noqa: E402
from thermowave.components.sink import Sink  # noqa: E402
from thermowave.components.source import Source  # noqa: E402
from thermowave.core.network import Network, NetworkState  # noqa: E402
from thermowave.fluids.real_fluid import CoolPropFluid  # noqa: E402

CO2 = CoolPropFluid(name="CO2")


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


def test_rejects_bad_PR():
    with pytest.raises(ValueError, match="PR"):
        SimpleHeater(name="h", T_out=500.0, PR=0.0)
    with pytest.raises(ValueError, match="PR"):
        SimpleHeater(name="h", T_out=500.0, PR=1.5)


def test_rejects_nonpositive_T_out():
    with pytest.raises(ValueError, match="T_out"):
        SimpleHeater(name="h", T_out=0.0)


def test_ports():
    h = SimpleHeater(name="h", T_out=500.0)
    assert h.ports() == {"in": "h.in", "out": "h.out"}


def test_residuals_hand_calc_supercritical_co2():
    # 250 bar is far above CO2's ~73.8 bar critical pressure -- exactly the
    # regime SimpleEvaporator can't handle (no saturation curve up there).
    P = 250.0e5
    T_out = 600.0 + 273.15
    h = SimpleHeater(name="h", T_out=T_out, PR=0.97276)
    h_out_target = CO2.enthalpy_pt(h.PR * P, T_out)
    P_out_expected = h.PR * P
    state = _FakeState(
        fluid=CO2,
        node_values={"h.in": (P, 4.0e5), "h.out": (P_out_expected, h_out_target)},
        mdots={"h.in": 2.0, "h.out": 2.0},
    )
    residuals = h.residuals(state)
    assert residuals == pytest.approx([0.0, 0.0, 0.0], abs=1e-6)


def test_network_solve_matches_target_outlet_temperature():
    net = Network(fluid=CO2)
    P_in = 257.0e5
    src = Source(name="src", P=P_in, T=433.63 + 273.15, mdot=100.0)
    T_out = 600.0 + 273.15
    heater = SimpleHeater(name="heater", PR=250.0 / 257.0, T_out=T_out)
    sink = Sink(name="sink")

    for c in (src, heater, sink):
        net.add_component(c)
    net.connect(src, "out", heater, "in")
    net.connect(heater, "out", sink, "in")

    result = net.solve(progress=False)
    assert result.converged

    state = result.state()
    metrics = heater.report_metrics(state)
    assert metrics["T_out [K]"] == pytest.approx(T_out, rel=1e-8)
    assert metrics["PR [-]"] == pytest.approx(250.0 / 257.0, rel=1e-6)
    assert metrics["power [W]"] > 0.0  # heating: duty is positive
