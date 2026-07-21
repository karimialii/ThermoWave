import math

import pytest

pytest.importorskip("cantera")

from thermowave.components.combustor import Combustor  # noqa: E402
from thermowave.components.sink import Sink  # noqa: E402
from thermowave.components.source import Source  # noqa: E402
from thermowave.core.network import Network  # noqa: E402
from thermowave.fluids.ideal_gas import IdealGasFluid  # noqa: E402

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def test_ports_and_category():
    comb = Combustor(name="cc1", mdot_fuel=0.02)
    assert comb.ports() == {"in": "cc1.in", "out": "cc1.out"}
    assert comb.report_category() == "combustor"


def test_free_parameters_empty_when_fuel_given():
    comb = Combustor(name="cc1", mdot_fuel=0.02)
    assert comb.free_parameters() == {}


def test_free_parameters_includes_fuel_when_omitted():
    comb = Combustor(name="cc1", mdot_fuel=None)
    assert "mdot_fuel" in comb.free_parameters()


def test_equilibrium_t_out_rises_with_fuel_flow():
    comb = Combustor(name="cc1", mdot_fuel=0.02)
    T_low = comb._equilibrium_T_out(T_in=500.0, P_in=3.0e5, mdot_air=1.0, mdot_fuel=0.01)
    T_high = comb._equilibrium_T_out(T_in=500.0, P_in=3.0e5, mdot_air=1.0, mdot_fuel=0.03)
    assert T_high > T_low > 500.0


def test_end_to_end_network_solve_with_fixed_fuel():
    src = Source(name="src", P=300000.0, T=500.0, mdot=1.0)
    comb = Combustor(name="cc1", PR=0.96, mdot_fuel=0.02)
    snk = Sink(name="snk")

    network = Network(fluid=AIR)
    for component in (src, comb, snk):
        network.add_component(component)
    network.connect(src, "out", comb, "in")
    network.connect(comb, "out", snk, "in")

    result = network.solve(tol=1e-6, max_iter=100)
    assert math.isclose(result.node_mdot["cc1.out"], 1.02, rel_tol=1e-9)
    T_out = AIR.temperature_ph(result.node_P["cc1.out"], result.node_h["cc1.out"])
    assert T_out > 500.0


def test_product_composition_returns_major_species_summing_near_one():
    from thermowave.core.network import NetworkState

    comb = Combustor(name="cc1", mdot_fuel=0.02)
    state = NetworkState(
        fluid=AIR,
        node_P={"cc1.in": 300000.0},
        node_h={"cc1.in": AIR.enthalpy_pt(300000.0, 500.0)},
        node_mdot={"cc1.in": 1.0},
    )
    products = comb.product_composition(state)
    assert "CO2" in products
    assert "H2O" in products
    assert "N2" in products
    assert math.isclose(sum(products.values()), 1.0, abs_tol=1e-3)
    for fraction in products.values():
        assert fraction >= 0.0


def test_product_composition_co2_and_h2o_rise_with_fuel_flow():
    from thermowave.core.network import NetworkState

    comb = Combustor(name="cc1", mdot_fuel=0.02)

    def _products(mdot_fuel):
        state = NetworkState(
            fluid=AIR,
            node_P={"cc1.in": 300000.0},
            node_h={"cc1.in": AIR.enthalpy_pt(300000.0, 500.0)},
            node_mdot={"cc1.in": 1.0},
            params={"cc1.mdot_fuel": mdot_fuel},
        )
        comb.mdot_fuel = mdot_fuel
        return comb.product_composition(state)

    products_lean = _products(0.01)
    products_rich = _products(0.03)
    assert products_rich["CO2"] > products_lean["CO2"]
    assert products_rich["H2O"] > products_lean["H2O"]


def test_report_metrics_includes_co2_and_h2o_mole_fractions():
    from thermowave.core.network import NetworkState

    comb = Combustor(name="cc1", mdot_fuel=0.02)
    P_in, T_in = 300000.0, 500.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    T_out = comb._equilibrium_T_out(T_in=T_in, P_in=P_in, mdot_air=1.0, mdot_fuel=0.02)
    P_out = comb.PR * P_in
    h_out = AIR.enthalpy_pt(P_out, T_out)
    state = NetworkState(
        fluid=AIR,
        node_P={"cc1.in": P_in, "cc1.out": P_out},
        node_h={"cc1.in": h_in, "cc1.out": h_out},
        node_mdot={"cc1.in": 1.0, "cc1.out": 1.02},
    )
    metrics = comb.report_metrics(state)
    assert "X_CO2 [-]" in metrics
    assert "X_H2O [-]" in metrics
    assert 0.0 < metrics["X_CO2 [-]"] < 1.0
    assert 0.0 < metrics["X_H2O [-]"] < 1.0


def test_end_to_end_network_solve_with_cantera_fluid():
    from thermowave.fluids.cantera_fluid import CanteraFluid

    cantera_air = CanteraFluid(name="cantera_air", composition="O2:0.21, N2:0.79")

    src = Source(name="src", P=300000.0, T=500.0, mdot=1.0)
    comb = Combustor(name="cc1", PR=0.96, mdot_fuel=0.02, fuel="CH4")
    snk = Sink(name="snk")

    network = Network(fluid=cantera_air)
    for component in (src, comb, snk):
        network.add_component(component)
    network.connect(src, "out", comb, "in")
    network.connect(comb, "out", snk, "in")

    result = network.solve(tol=1e-6, max_iter=200)
    assert math.isclose(result.node_mdot["cc1.out"], 1.02, rel_tol=1e-9)

    # Composition propagates downstream of the combustor: the outlet node's
    # fluid is the reacted product mixture (CO2/H2O present), not the plain
    # pre-combustion air the network was built with.
    outlet_fluid = result.node_fluid["cc1.out"]
    assert outlet_fluid is not cantera_air
    T_out = outlet_fluid.temperature_ph(result.node_P["cc1.out"], result.node_h["cc1.out"])
    assert T_out > 500.0

    # Reading the same (P, h) through the *wrong* (pre-combustion) air
    # composition gives a materially different answer -- proof the outlet
    # node's fluid is actually a different composition, not just a
    # differently-named alias for the same one.
    T_out_via_air = cantera_air.temperature_ph(
        result.node_P["cc1.out"], result.node_h["cc1.out"]
    )
    assert not math.isclose(T_out, T_out_via_air, rel_tol=1e-3)


class _FakeHeatPath:
    def __init__(self, Q_value):
        self._Q_value = Q_value

    def Q(self, state):
        return self._Q_value


def test_combustor_heat_path_defaults_to_none():
    comb = Combustor(name="cc1", mdot_fuel=0.02)
    assert comb.heat_path is None


def test_combustor_heat_path_adds_q_over_outlet_mdot_to_energy_residual():
    from thermowave.core.network import NetworkState

    mdot_in, mdot_fuel = 1.0, 0.02
    mdot_out = mdot_in + mdot_fuel
    P_in, T_in = 300000.0, 500.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    P_out, h_out = 0.96 * P_in, h_in  # arbitrary outlet guess, shared by both cases

    state = NetworkState(
        fluid=AIR,
        node_P={"cc1.in": P_in, "cc1.out": P_out},
        node_h={"cc1.in": h_in, "cc1.out": h_out},
        node_mdot={"cc1.in": mdot_in, "cc1.out": mdot_out},
    )

    baseline = Combustor(name="cc1", mdot_fuel=mdot_fuel)
    _, baseline_energy_residual, _ = baseline.residuals(state)

    Q_loss = 800.0
    lossy = Combustor(name="cc1", mdot_fuel=mdot_fuel)
    lossy.heat_path = _FakeHeatPath(Q_value=Q_loss)
    _, lossy_energy_residual, _ = lossy.residuals(state)

    assert math.isclose(
        lossy_energy_residual - baseline_energy_residual, Q_loss / mdot_out, rel_tol=1e-9
    )


def test_combustor_report_metrics_includes_q_loss():
    from thermowave.core.network import NetworkState

    comb = Combustor(name="cc1", mdot_fuel=0.02)
    P_in, T_in = 300000.0, 500.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    state = NetworkState(
        fluid=AIR,
        node_P={"cc1.in": P_in, "cc1.out": 0.96 * P_in},
        node_h={"cc1.in": h_in, "cc1.out": h_in},
        node_mdot={"cc1.in": 1.0, "cc1.out": 1.02},
    )
    assert comb.report_metrics(state)["Q_loss [W]"] == 0.0

    comb.heat_path = _FakeHeatPath(Q_value=456.0)
    assert comb.report_metrics(state)["Q_loss [W]"] == 456.0


def test_import_error_message_when_cantera_missing(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "cantera" or name.startswith("cantera"):
            raise ImportError("simulated missing cantera")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="cantera"):
        Combustor(name="cc1", mdot_fuel=0.02)
