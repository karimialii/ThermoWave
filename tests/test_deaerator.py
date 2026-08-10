import math

import pytest

pytest.importorskip("CoolProp")

from thermowave.components.deaerator import Deaerator  # noqa: E402
from thermowave.components.sink import Sink  # noqa: E402
from thermowave.components.source import Source  # noqa: E402
from thermowave.core.network import Network  # noqa: E402
from thermowave.fluids.ideal_gas import IdealGasFluid  # noqa: E402
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


def test_rejects_fewer_than_two_inlets():
    with pytest.raises(ValueError, match="n_inlets"):
        Deaerator(name="fwt", n_inlets=1)


def test_rejects_non_two_phase_fluid():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    d = Deaerator(name="fwt", n_inlets=3)
    state = _FakeState(
        fluid=air,
        node_values={
            "fwt.in0": (5.0e5, 1.0e5), "fwt.in1": (5.0e5, 1.0e5),
            "fwt.in2": (5.0e5, 1.0e5), "fwt.out": (5.0e5, 1.0e5),
        },
        mdots={"fwt.in0": 1.0, "fwt.in1": 1.0, "fwt.in2": 1.0, "fwt.out": 3.0},
    )
    with pytest.raises(ValueError, match="two-phase-capable"):
        d.residuals(state)


def test_ports_has_n_inlets_plus_one_outlet():
    d = Deaerator(name="fwt", n_inlets=3)
    assert set(d.ports()) == {"in0", "in1", "in2", "out"}


def test_residuals_zero_at_saturated_liquid_outlet_hand_calc():
    d = Deaerator(name="fwt", n_inlets=3)
    P = 5.0e5
    h_f = WATER.saturated_liquid_enthalpy(P)
    # Inlet enthalpies don't matter to the residuals themselves (the outlet
    # is pinned directly to h_f, not to a mass-weighted average of these) --
    # only mass conservation and the saturation pin are checked.
    state = _FakeState(
        fluid=WATER,
        node_values={
            "fwt.in0": (P, WATER.enthalpy_pt(P, 400.0)),
            "fwt.in1": (P, WATER.enthalpy_pt(P, 420.0)),
            "fwt.in2": (P, h_f),
            "fwt.out": (P, h_f),
        },
        mdots={"fwt.in0": 1.0, "fwt.in1": 2.0, "fwt.in2": 0.5, "fwt.out": 3.5},
    )
    residuals = d.residuals(state)
    assert len(residuals) == 3
    for r in residuals:
        assert math.isclose(r, 0.0, abs_tol=1e-6)


def test_residuals_nonzero_when_outlet_off_saturation_or_mass_unbalanced():
    d = Deaerator(name="fwt", n_inlets=2)
    P = 5.0e5
    h_f = WATER.saturated_liquid_enthalpy(P)
    state = _FakeState(
        fluid=WATER,
        node_values={
            "fwt.in0": (P, WATER.enthalpy_pt(P, 400.0)),
            "fwt.in1": (P, WATER.enthalpy_pt(P, 420.0)),
            "fwt.out": (P, h_f + 5000.0),  # deliberately superheated, not saturated liquid
        },
        mdots={"fwt.in0": 1.0, "fwt.in1": 1.0, "fwt.out": 1.5},  # deliberately mass-unbalanced
    )
    momentum, energy, mass = d.residuals(state)
    assert math.isclose(momentum, 0.0, abs_tol=1e-6)  # P_out matches P_ref (in0)
    assert not math.isclose(energy, 0.0, abs_tol=1.0)
    assert not math.isclose(mass, 0.0, abs_tol=1e-6)


def test_report_metrics_reports_saturation_temperature():
    d = Deaerator(name="fwt", n_inlets=2)
    P = 5.0e5
    h_f = WATER.saturated_liquid_enthalpy(P)
    state = _FakeState(
        fluid=WATER,
        node_values={
            "fwt.in0": (P, WATER.enthalpy_pt(P, 400.0)),
            "fwt.in1": (P, WATER.enthalpy_pt(P, 420.0)),
            "fwt.out": (P, h_f),
        },
        mdots={"fwt.in0": 1.0, "fwt.in1": 2.0, "fwt.out": 3.0},
    )
    metrics = d.report_metrics(state)
    assert math.isclose(metrics["T_sat [K]"], WATER.saturation_temperature(P), rel_tol=1e-9)
    assert math.isclose(metrics["mdot_out [kg/s]"], 3.0, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# End-to-end: a real Network.solve() with Deaerator mixing two subcooled
# streams and settling at saturated liquid -- proves the equilibrium
# residual actually closes the system, not just that its shape is right.
# ---------------------------------------------------------------------------


def test_deaerator_end_to_end_solve_settles_at_saturated_liquid():
    P = 5.0e5
    src0 = Source(name="src0", P=P, T=400.0, mdot=1.0)
    src1 = Source(name="src1", P=P, T=420.0, mdot=2.0)
    d = Deaerator(name="fwt", n_inlets=2)
    snk = Sink(name="snk")

    network = Network(fluid=WATER)
    for c in (src0, src1, d, snk):
        network.add_component(c)
    network.connect(src0, "out", d, "in0")
    network.connect(src1, "out", d, "in1")
    network.connect(d, "out", snk, "in")

    result = network.solve(tol=1e-9, max_iter=100)
    assert result.converged

    P_out, h_out = result.node_P["fwt.out"], result.node_h["fwt.out"]
    h_f = WATER.saturated_liquid_enthalpy(P_out)
    assert math.isclose(h_out, h_f, rel_tol=1e-6)
    assert math.isclose(result.node_mdot["fwt.out"], 3.0, rel_tol=1e-9)
