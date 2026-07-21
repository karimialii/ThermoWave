import math

from thermowave.components.compressor import Compressor
from thermowave.components.heat_transfer import Conduction, Convection, ThermalMass
from thermowave.components.pipe import Pipe
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)
GAMMA = 1005.0 / (1005.0 - 287.05)
T_AMBIENT = 288.15


def _build_network_with_casings(N_shaft=65000.0):
    """Source -> Compressor -> Pipe(heater) -> Turbine -> Sink, the same
    fixed-N recipe examples/19 uses, with each machine's fluid convecting
    heat into its own casing ThermalMass, the two casings conducting through
    a shaft ThermalMass between them, and both casings convecting to
    ambient -- "turbine loses heat to ambient through its casing, and heat
    conducts through the turbine case and shaft to the compressor case."
    """
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    comp = Compressor(name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=N_shaft)
    heater = Pipe(name="heater", L=1.0, D=0.1, f=0.0, n_elem=1, heat_loss=-431000.0)
    turb = Turbine(name="turb", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=N_shaft)
    snk = Sink(name="snk")

    comp_casing = ThermalMass(name="comp_casing", thermal_capacitance=200.0, T0=300.0)
    turb_casing = ThermalMass(name="turb_casing", thermal_capacitance=200.0, T0=300.0)
    shaft_mass = ThermalMass(name="shaft_mass", thermal_capacitance=50.0, T0=300.0)

    conv_comp_fluid = Convection(
        name="conv_comp_fluid", a=(comp, "out"), b=comp_casing, h=50.0, A=0.3
    )
    conv_turb_fluid = Convection(
        name="conv_turb_fluid", a=(turb, "out"), b=turb_casing, h=50.0, A=0.3
    )
    cond_turb_shaft = Conduction(
        name="cond_turb_shaft", a=turb_casing, b=shaft_mass, k=15.0, A=0.01, L=0.5
    )
    cond_shaft_comp = Conduction(
        name="cond_shaft_comp", a=shaft_mass, b=comp_casing, k=15.0, A=0.01, L=0.5
    )
    conv_turb_ambient = Convection(
        name="conv_turb_ambient", a=turb_casing, b=T_AMBIENT, h=10.0, A=1.0
    )
    conv_comp_ambient = Convection(
        name="conv_comp_ambient", a=comp_casing, b=T_AMBIENT, h=10.0, A=1.0
    )

    comp.heat_path = conv_comp_fluid
    turb.heat_path = conv_turb_fluid
    comp_casing.heat_sources = [
        (conv_comp_fluid, 1.0), (cond_shaft_comp, 1.0), (conv_comp_ambient, -1.0),
    ]
    turb_casing.heat_sources = [
        (conv_turb_fluid, 1.0), (cond_turb_shaft, -1.0), (conv_turb_ambient, -1.0),
    ]
    shaft_mass.heat_sources = [(cond_turb_shaft, 1.0), (cond_shaft_comp, -1.0)]

    network = Network(fluid=AIR)
    for component in (
        src, comp, heater, turb, snk,
        comp_casing, turb_casing, shaft_mass,
        conv_comp_fluid, conv_turb_fluid, cond_turb_shaft, cond_shaft_comp,
        conv_turb_ambient, conv_comp_ambient,
    ):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", heater, "in")
    network.connect(heater, "out", turb, "in")
    network.connect(turb, "out", snk, "in")
    return network, comp, turb, comp_casing, turb_casing, shaft_mass


def test_casing_temperatures_rise_monotonically_from_cold_start():
    network, comp, turb, comp_casing, turb_casing, shaft_mass = _build_network_with_casings()

    # solve_transient()'s default initial condition is a genuine steady
    # state (every differential state's own derivative already zero at
    # t=0), which for these casings means already at thermal equilibrium --
    # not a cold start. Build a real t=0 SolveResult the normal way, then
    # override just the three thermal-mass states to a uniform cold value
    # before handing it to solve_transient() as `initial`, so this test
    # actually exercises the casings warming up over time.
    cold_initial = network.solve(tol=1e-8, max_iter=400, damping=0.3)
    for mass_name in ("comp_casing", "turb_casing", "shaft_mass"):
        cold_initial.params[f"{mass_name}.T"] = 300.0

    history = network.solve_transient(
        duration=200.0, dt=10.0, tol=1e-8, max_iter=400, damping=0.3, initial=cold_initial,
    )

    turb_T = history.diff_history["turb_casing.T"]
    comp_T = history.diff_history["comp_casing.T"]
    shaft_T = history.diff_history["shaft_mass.T"]

    for series in (turb_T, comp_T, shaft_T):
        assert series[0] == 300.0
        assert all(b >= a - 1e-9 for a, b in zip(series, series[1:])), series
        assert series[-1] > series[0]

    # The turbine casing (hot-fluid side) should end up hotter than the
    # compressor casing (cooler-fluid side) once heat has propagated.
    assert turb_T[-1] > comp_T[-1]


def test_turbine_outlet_temperature_is_lower_with_heat_loss_than_without():
    network, comp, turb, comp_casing, turb_casing, shaft_mass = _build_network_with_casings()
    history = network.solve_transient(
        duration=200.0, dt=10.0, tol=1e-8, max_iter=400, damping=0.3,
    )
    last_step = history.steps[-1]
    T_out_with_loss = AIR.temperature_ph(
        last_step.node_P["turb.out"], last_step.node_h["turb.out"]
    )

    baseline_network, baseline_comp, baseline_turb, *_ = _build_network_with_casings()
    # No heat paths at all on the baseline -- fully adiabatic turbomachinery,
    # exactly like every other example in this codebase.
    baseline_turb.heat_path = None
    baseline_comp.heat_path = None
    baseline_result = baseline_network.solve(tol=1e-8, max_iter=400, damping=0.3)
    T_out_baseline = AIR.temperature_ph(
        baseline_result.node_P["turb.out"], baseline_result.node_h["turb.out"]
    )

    assert T_out_with_loss < T_out_baseline


def test_thermal_masses_reach_zero_net_heat_at_network_solve_steady_state():
    network, comp, turb, comp_casing, turb_casing, shaft_mass = _build_network_with_casings()
    result = network.solve(tol=1e-9, max_iter=400, damping=0.3)

    from thermowave.core.network import NetworkState

    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params, node_fluid=result.node_fluid,
    )
    for mass in (comp_casing, turb_casing, shaft_mass):
        Q_net = mass._net_heat(state)
        assert math.isclose(Q_net, 0.0, abs_tol=1e-4)
