import math

from thermowave.components.pipe import Pipe
from thermowave.components.simple_combustor import SimpleCombustor
from thermowave.core.network import Network
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _FakeState:
    def __init__(self, fluid, mdot, node_values, params=None):
        self.fluid = fluid
        self._mdot = mdot
        self._node_values = node_values
        self.params = params or {}

    def node(self, name):
        return self._node_values[name]

    def fluid_at(self, name):
        return self.fluid

    def mdot(self, name):
        return self._mdot[name]

    def param(self, name):
        return self.params[name]


def test_ports_and_category():
    comb = SimpleCombustor(name="cc1", LHV=50e6, mdot_fuel=0.02)
    assert comb.ports() == {"in": "cc1.in", "out": "cc1.out"}
    assert comb.report_category() == "combustor"


def test_free_parameters_empty_when_fuel_given():
    comb = SimpleCombustor(name="cc1", LHV=50e6, mdot_fuel=0.02)
    assert comb.free_parameters() == {}


def test_free_parameters_includes_fuel_when_omitted():
    comb = SimpleCombustor(name="cc1", LHV=50e6, mdot_fuel=None)
    assert "mdot_fuel" in comb.free_parameters()


def test_residuals_match_hand_calc_with_fixed_fuel():
    P_in, T_in = 300000.0, 500.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    mdot_in = 1.0
    mdot_fuel = 0.02
    LHV, eta, PR = 50e6, 0.98, 0.96

    comb = SimpleCombustor(name="cc1", LHV=LHV, PR=PR, efficiency=eta, mdot_fuel=mdot_fuel)

    Q = mdot_fuel * LHV * eta
    mdot_out = mdot_in + mdot_fuel
    h_out = (mdot_in * h_in + Q) / mdot_out
    P_out = PR * P_in

    state = _FakeState(
        fluid=AIR,
        mdot={"cc1.in": mdot_in, "cc1.out": mdot_out},
        node_values={"cc1.in": (P_in, h_in), "cc1.out": (P_out, h_out)},
    )
    residuals = comb.residuals(state)
    assert len(residuals) == 3
    for r in residuals:
        assert math.isclose(r, 0.0, abs_tol=1e-3)


def test_residuals_use_free_fuel_from_state_param_when_omitted():
    P_in, T_in = 300000.0, 500.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    mdot_in = 1.0
    mdot_fuel = 0.03
    LHV, PR = 45e6, 0.95

    comb = SimpleCombustor(name="cc1", LHV=LHV, PR=PR, mdot_fuel=None)

    Q = mdot_fuel * LHV * 1.0
    mdot_out = mdot_in + mdot_fuel
    h_out = (mdot_in * h_in + Q) / mdot_out
    P_out = PR * P_in

    state = _FakeState(
        fluid=AIR,
        mdot={"cc1.in": mdot_in, "cc1.out": mdot_out},
        node_values={"cc1.in": (P_in, h_in), "cc1.out": (P_out, h_out)},
        params={"cc1.mdot_fuel": mdot_fuel},
    )
    residuals = comb.residuals(state)
    for r in residuals:
        assert math.isclose(r, 0.0, abs_tol=1e-3)


def test_report_metrics_includes_fuel_pr_and_t_out():
    P_in, T_in = 300000.0, 500.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    mdot_in = 1.0
    mdot_fuel = 0.02
    LHV, PR = 50e6, 0.96
    comb = SimpleCombustor(name="cc1", LHV=LHV, PR=PR, mdot_fuel=mdot_fuel)

    Q = mdot_fuel * LHV
    mdot_out = mdot_in + mdot_fuel
    h_out = (mdot_in * h_in + Q) / mdot_out
    P_out = PR * P_in
    T_out_expected = AIR.temperature_ph(P_out, h_out)

    state = _FakeState(
        fluid=AIR,
        mdot={"cc1.in": mdot_in, "cc1.out": mdot_out},
        node_values={"cc1.in": (P_in, h_in), "cc1.out": (P_out, h_out)},
    )
    metrics = comb.report_metrics(state)
    assert math.isclose(metrics["mdot_fuel [kg/s]"], mdot_fuel)
    assert math.isclose(metrics["PR [-]"], PR)
    assert math.isclose(metrics["T_out [K]"], T_out_expected, rel_tol=1e-9)
    assert math.isclose(metrics["power [W]"], Q, rel_tol=1e-6)


def test_end_to_end_network_solve_with_fixed_fuel():
    src = Source(name="src", P=300000.0, T=500.0, mdot=1.0)
    comb = SimpleCombustor(name="cc1", LHV=50e6, PR=0.96, mdot_fuel=0.02)
    snk = Sink(name="snk")

    network = Network(fluid=AIR)
    for component in (src, comb, snk):
        network.add_component(component)
    network.connect(src, "out", comb, "in")
    network.connect(comb, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=100)
    assert math.isclose(result.node_mdot["cc1.out"], 1.02, rel_tol=1e-9)
    T_out = AIR.temperature_ph(result.node_P["cc1.out"], result.node_h["cc1.out"])
    assert T_out > 500.0


# --- use_fuel_port=True (real two-inlet mode) -----------------------------


def test_ports_includes_fuel_in_when_use_fuel_port():
    comb = SimpleCombustor(name="cc1", LHV=50e6, use_fuel_port=True)
    assert comb.ports() == {
        "in": "cc1.in",
        "out": "cc1.out",
        "fuel_in": "cc1.fuel_in",
    }


def test_free_parameters_empty_when_use_fuel_port():
    comb = SimpleCombustor(name="cc1", LHV=50e6, use_fuel_port=True)
    assert comb.free_parameters() == {}


def test_residuals_use_fuel_port_mdot_and_include_fuel_enthalpy():
    P_in, T_in = 300000.0, 500.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    P_fuel, T_fuel = 350000.0, 300.0
    h_fuel = AIR.enthalpy_pt(P_fuel, T_fuel)
    mdot_in = 1.0
    mdot_fuel = 0.02
    LHV, eta, PR = 50e6, 0.98, 0.96

    comb = SimpleCombustor(name="cc1", LHV=LHV, PR=PR, efficiency=eta, use_fuel_port=True)

    Q = mdot_fuel * LHV * eta
    mdot_out = mdot_in + mdot_fuel
    h_out = (mdot_in * h_in + mdot_fuel * h_fuel + Q) / mdot_out
    P_out = PR * P_in

    state = _FakeState(
        fluid=AIR,
        mdot={"cc1.in": mdot_in, "cc1.out": mdot_out, "cc1.fuel_in": mdot_fuel},
        node_values={
            "cc1.in": (P_in, h_in),
            "cc1.out": (P_out, h_out),
            "cc1.fuel_in": (P_fuel, h_fuel),
        },
    )
    residuals = comb.residuals(state)
    assert len(residuals) == 3
    for r in residuals:
        assert math.isclose(r, 0.0, abs_tol=1e-3)


def test_end_to_end_network_solve_with_fuel_port_and_fixed_fuel_source():
    src = Source(name="src", P=300000.0, T=500.0, mdot=1.0)
    fuel_src = Source(name="fuel_src", P=350000.0, T=300.0, mdot=0.02)
    fuel_pipe = Pipe(name="fuel_pipe", L=1.0, D=0.01, f=0.02, n_elem=1)
    comb = SimpleCombustor(name="cc1", LHV=50e6, PR=0.96, use_fuel_port=True)
    snk = Sink(name="snk")

    network = Network(fluid=AIR)
    for component in (src, fuel_src, fuel_pipe, comb, snk):
        network.add_component(component)
    network.connect(src, "out", comb, "in")
    network.connect(fuel_src, "out", fuel_pipe, "in")
    network.connect(fuel_pipe, "out", comb, "fuel_in")
    network.connect(comb, "out", snk, "in")

    result = network.solve(tol=1e-8, max_iter=100)
    assert math.isclose(result.node_mdot["cc1.out"], 1.02, rel_tol=1e-6)
    T_out = AIR.temperature_ph(result.node_P["cc1.out"], result.node_h["cc1.out"])
    assert T_out > 500.0


class _FakeHeatPath:
    def __init__(self, Q_value):
        self._Q_value = Q_value

    def Q(self, state):
        return self._Q_value


def test_simple_combustor_heat_path_defaults_to_none():
    comb = SimpleCombustor(name="cc1", LHV=50e6, mdot_fuel=0.02)
    assert comb.heat_path is None


def test_simple_combustor_heat_path_subtracts_q_from_power_valued_energy_residual():
    P_in, T_in = 300000.0, 500.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    mdot_in = 1.0
    mdot_fuel = 0.02
    LHV, eta, PR = 50e6, 0.98, 0.96

    Q = mdot_fuel * LHV * eta
    mdot_out = mdot_in + mdot_fuel
    h_out = (mdot_in * h_in + Q) / mdot_out
    P_out = PR * P_in

    state = _FakeState(
        fluid=AIR,
        mdot={"cc1.in": mdot_in, "cc1.out": mdot_out},
        node_values={"cc1.in": (P_in, h_in), "cc1.out": (P_out, h_out)},
    )

    baseline = SimpleCombustor(name="cc1", LHV=LHV, PR=PR, efficiency=eta, mdot_fuel=mdot_fuel)
    _, baseline_energy_residual, _ = baseline.residuals(state)

    Q_loss = 500.0
    lossy = SimpleCombustor(name="cc1", LHV=LHV, PR=PR, efficiency=eta, mdot_fuel=mdot_fuel)
    lossy.heat_path = _FakeHeatPath(Q_value=Q_loss)
    _, lossy_energy_residual, _ = lossy.residuals(state)

    # energy_residual = mdot_out*h_out - (... + Q - Q_loss): increasing
    # Q_loss increases the residual (less energy is actually available to
    # heat the fluid, so the same h_out guess now overshoots more).
    assert math.isclose(
        lossy_energy_residual - baseline_energy_residual, Q_loss, rel_tol=1e-9
    )


def test_simple_combustor_report_metrics_includes_q_loss():
    P_in, T_in = 300000.0, 500.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    mdot_in, mdot_fuel = 1.0, 0.02
    state = _FakeState(
        fluid=AIR,
        mdot={"cc1.in": mdot_in, "cc1.out": mdot_in + mdot_fuel},
        node_values={"cc1.in": (P_in, h_in), "cc1.out": (0.96 * P_in, h_in)},
    )
    comb = SimpleCombustor(name="cc1", LHV=50e6, mdot_fuel=mdot_fuel)
    assert comb.report_metrics(state)["Q_loss [W]"] == 0.0

    comb.heat_path = _FakeHeatPath(Q_value=321.0)
    assert comb.report_metrics(state)["Q_loss [W]"] == 321.0


def test_simple_combustor_rejects_pr_out_of_range():
    try:
        SimpleCombustor(name="cc1", LHV=45.0e6, PR=1.2)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_simple_combustor_rejects_non_positive_efficiency():
    try:
        SimpleCombustor(name="cc1", LHV=45.0e6, efficiency=0.0)
        assert False, "expected ValueError"
    except ValueError:
        pass
