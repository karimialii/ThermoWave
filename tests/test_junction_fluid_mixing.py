import math

import pytest

from thermowave.components.junction import Junction
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)
OTHER_IDEAL = IdealGasFluid(name="other", R=300.0, cp=1100.0)


class _FakeState:
    def __init__(self, mdot_values):
        self._mdot_values = mdot_values

    def mdot(self, name: str) -> float:
        return self._mdot_values[name]


def test_junction_fluid_flow_pairs_is_empty():
    # Junction resolves outlet fluids entirely through merge_fluids(),
    # never the generic single-inlet-per-pair path -- see its own
    # fluid_flow_pairs() docstring.
    j = Junction(name="j1", n_inlets=2, n_outlets=2)
    assert j.fluid_flow_pairs() == []


def test_merge_fluids_returns_none_when_not_every_inlet_is_resolved_yet():
    j = Junction(name="j1", n_inlets=2, n_outlets=1)
    state = _FakeState(mdot_values={})
    assert j.merge_fluids(state, {"in0": AIR}) is None


def test_merge_fluids_passes_through_when_every_inlet_is_the_same_fluid():
    j = Junction(name="j1", n_inlets=2, n_outlets=2)
    state = _FakeState(mdot_values={"j1.in0": 1.0, "j1.in1": 1.0})
    merged = j.merge_fluids(state, {"in0": AIR, "in1": AIR})
    assert merged == {"out0": AIR, "out1": AIR}


def test_merge_fluids_falls_back_to_first_inlet_when_not_blendable():
    # Plain IdealGasFluid exposes neither mass_fractions() nor mechanism --
    # no principled composition blend is possible, so this falls back to
    # the documented "first inlet wins" simplification.
    j = Junction(name="j1", n_inlets=2, n_outlets=1)
    state = _FakeState(mdot_values={"j1.in0": 1.0, "j1.in1": 1.0})
    merged = j.merge_fluids(state, {"in0": AIR, "in1": OTHER_IDEAL})
    assert merged == {"out0": AIR}


cantera = pytest.importorskip("cantera")


def test_merge_fluids_falls_back_to_first_inlet_for_different_mechanisms():
    from thermowave.fluids.cantera_fluid import CanteraFluid

    air1 = CanteraFluid(name="air1", composition="O2:0.21, N2:0.79", mechanism="gri30.yaml")
    air2 = CanteraFluid(name="air2", composition="O2:0.21, N2:0.79", mechanism="h2o2.yaml")
    j = Junction(name="j1", n_inlets=2, n_outlets=1)
    state = _FakeState(mdot_values={"j1.in0": 1.0, "j1.in1": 1.0})
    merged = j.merge_fluids(state, {"in0": air1, "in1": air2})
    assert merged == {"out0": air1}


def test_merge_fluids_blends_two_cantera_compositions_by_mass():
    from thermowave.fluids.cantera_fluid import CanteraFluid

    fluid_a = CanteraFluid(name="a", composition="N2:1.0", mechanism="gri30.yaml")
    fluid_b = CanteraFluid(name="b", composition="O2:1.0", mechanism="gri30.yaml")
    j = Junction(name="j1", n_inlets=2, n_outlets=1)
    mdot_a, mdot_b = 3.0, 1.0
    state = _FakeState(mdot_values={"j1.in0": mdot_a, "j1.in1": mdot_b})

    merged = j.merge_fluids(state, {"in0": fluid_a, "in1": fluid_b})
    mixed = merged["out0"]
    assert mixed is not fluid_a
    assert mixed is not fluid_b

    Y = mixed.mass_fractions()
    expected_N2 = mdot_a / (mdot_a + mdot_b)
    expected_O2 = mdot_b / (mdot_a + mdot_b)
    assert math.isclose(Y.get("N2", 0.0), expected_N2, rel_tol=1e-9)
    assert math.isclose(Y.get("O2", 0.0), expected_O2, rel_tol=1e-9)


def test_merge_fluids_blend_round_trips_enthalpy_temperature():
    from thermowave.fluids.cantera_fluid import CanteraFluid

    fluid_a = CanteraFluid(name="a", composition="N2:1.0", mechanism="gri30.yaml")
    fluid_b = CanteraFluid(name="b", composition="O2:1.0", mechanism="gri30.yaml")
    j = Junction(name="j1", n_inlets=2, n_outlets=1)
    state = _FakeState(mdot_values={"j1.in0": 1.0, "j1.in1": 1.0})
    mixed = j.merge_fluids(state, {"in0": fluid_a, "in1": fluid_b})["out0"]

    P, T_probe = 200000.0, 600.0
    h = mixed.enthalpy_pt(P, T_probe)
    T_recovered = mixed.temperature_ph(P, h)
    assert math.isclose(T_recovered, T_probe, rel_tol=1e-6)


def test_merge_fluids_all_zero_inflow_does_not_raise():
    from thermowave.fluids.cantera_fluid import CanteraFluid

    fluid_a = CanteraFluid(name="a", composition="N2:1.0", mechanism="gri30.yaml")
    fluid_b = CanteraFluid(name="b", composition="O2:1.0", mechanism="gri30.yaml")
    j = Junction(name="j1", n_inlets=2, n_outlets=1)
    state = _FakeState(mdot_values={"j1.in0": 0.0, "j1.in1": 0.0})
    merged = j.merge_fluids(state, {"in0": fluid_a, "in1": fluid_b})
    assert merged["out0"] is fluid_a


def test_merge_fluids_blends_combustor_product_fluids_not_just_cantera_fluid():
    # Regression guard: Combustor's own outlet_fluid() returns a
    # _CanteraCompositionFluid, not a CanteraFluid instance -- the gate must
    # be duck-typed (mass_fractions() + mechanism), not isinstance(...,
    # CanteraFluid), or real mixing would never trigger for the realistic
    # "two combustor exhaust streams merge" case this feature exists for.
    from thermowave.components.combustor import Combustor
    from thermowave.core.network import NetworkState
    from thermowave.fluids.cantera_fluid import CanteraFluid

    cantera_air = CanteraFluid(name="cantera_air", composition="O2:0.21, N2:0.79")
    comb_a = Combustor(name="cc_a", mdot_fuel=0.04, fuel="CH4")
    comb_b = Combustor(name="cc_b", mdot_fuel=0.01, fuel="CH4")
    P_in, T_in = 300000.0, 500.0
    h_in = cantera_air.enthalpy_pt(P_in, T_in)
    state = NetworkState(
        fluid=cantera_air,
        node_P={"cc_a.in": P_in, "cc_b.in": P_in},
        node_h={"cc_a.in": h_in, "cc_b.in": h_in},
        node_mdot={"cc_a.in": 1.0, "cc_b.in": 1.0},
    )
    product_a = comb_a.outlet_fluid(state, ("in", "out"), cantera_air)
    product_b = comb_b.outlet_fluid(state, ("in", "out"), cantera_air)
    assert type(product_a).__name__ == "_CanteraCompositionFluid"

    j = Junction(name="j1", n_inlets=2, n_outlets=1)
    jstate = _FakeState(mdot_values={"j1.in0": 1.0, "j1.in1": 1.0})
    merged = j.merge_fluids(jstate, {"in0": product_a, "in1": product_b})
    mixed = merged["out0"]
    assert mixed is not product_a
    assert mixed is not product_b
    Y_mixed = mixed.mass_fractions()
    Y_a, Y_b = product_a.mass_fractions(), product_b.mass_fractions()
    assert math.isclose(Y_mixed["CO2"], 0.5 * (Y_a["CO2"] + Y_b["CO2"]), rel_tol=1e-9)


# ---------------------------------------------------------------------------
# End-to-end: two combustor exhaust streams (genuinely different products)
# merging at a Junction
# ---------------------------------------------------------------------------


def test_junction_end_to_end_mixes_two_combustor_product_streams():
    from thermowave.components.combustor import Combustor
    from thermowave.components.sink import Sink
    from thermowave.components.source import Source
    from thermowave.core.network import Network
    from thermowave.fluids.cantera_fluid import CanteraFluid

    cantera_air = CanteraFluid(name="cantera_air", composition="O2:0.21, N2:0.79")

    src_a = Source(name="src_a", P=300000.0, T=500.0, mdot=1.0)
    comb_a = Combustor(name="cc_a", PR=0.96, mdot_fuel=0.04, fuel="CH4")  # richer
    src_b = Source(name="src_b", P=300000.0, T=500.0, mdot=1.0)
    comb_b = Combustor(name="cc_b", PR=0.96, mdot_fuel=0.01, fuel="CH4")  # leaner
    junction = Junction(name="j1", n_inlets=2, n_outlets=1)
    snk = Sink(name="snk")

    network = Network(fluid=cantera_air)
    for c in (src_a, comb_a, src_b, comb_b, junction, snk):
        network.add_component(c)
    network.connect(src_a, "out", comb_a, "in")
    network.connect(comb_a, "out", junction, "in0")
    network.connect(src_b, "out", comb_b, "in")
    network.connect(comb_b, "out", junction, "in1")
    network.connect(junction, "out0", snk, "in")

    result = network.solve(tol=1e-6, max_iter=400)
    assert result.converged

    mixed_fluid = result.node_fluid["j1.out0"]
    assert mixed_fluid is not cantera_air

    Y_mixed = mixed_fluid.mass_fractions()
    mdot_a = result.node_mdot["cc_a.out"]
    mdot_b = result.node_mdot["cc_b.out"]

    fluid_a = result.node_fluid["cc_a.out"]
    fluid_b = result.node_fluid["cc_b.out"]
    Y_a = fluid_a.mass_fractions()
    Y_b = fluid_b.mass_fractions()

    all_species = set(Y_a) | set(Y_b)
    for species in all_species:
        expected = (mdot_a * Y_a.get(species, 0.0) + mdot_b * Y_b.get(species, 0.0)) / (
            mdot_a + mdot_b
        )
        assert math.isclose(Y_mixed.get(species, 0.0), expected, abs_tol=1e-9)

    # CO2 mass fraction from the richer combustor must exceed the leaner one,
    # and the mixed stream must land strictly between them.
    assert Y_a["CO2"] > Y_b["CO2"]
    assert Y_b["CO2"] < Y_mixed["CO2"] < Y_a["CO2"]
