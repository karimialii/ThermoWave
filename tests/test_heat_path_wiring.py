"""Network.add_heat_path() — the one-call replacement for hand-wiring a
Convection/Conduction/Radiation onto both of its endpoints.

The reference test is test_add_heat_path_reproduces_hand_wiring_exactly: the
same physical network built both ways must produce byte-identical signs and
the same converged solution. Everything else here pins down the parts that
were previously silent failures (wrong sign, forgotten registration, bad
port name).
"""

import math

import pytest

from thermowave.components.compressor import Compressor
from thermowave.components.heat_transfer import (
    Conduction,
    Convection,
    Radiation,
    ThermalMass,
    heat_loss_watts,
    normalized_heat_paths,
)
from thermowave.components.pipe import Pipe
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.exceptions import NetworkTopologyError
from thermowave.core.network import Network
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)
GAMMA = 1005.0 / (1005.0 - 287.05)
T_AMBIENT = 288.15
_COMP_MAP = "tests/fixtures/simple_compressor_map.cop"
_TURB_MAP = "tests/fixtures/simple_turbine_map.tur"


def _flow_parts(N_shaft=65000.0):
    """Source -> Compressor -> Pipe(heater) -> Turbine -> Sink, the same
    fixed-N recipe tests/test_heat_transfer_integration.py uses."""
    return (
        Source(name="src", P=101325.0, T=288.15, mdot=0.63),
        Compressor(name="comp", map_path=_COMP_MAP, gamma=GAMMA, N=N_shaft),
        Pipe(name="heater", L=1.0, D=0.1, f=0.0, n_elem=1, heat_loss=-431000.0),
        Turbine(name="turb", map_path=_TURB_MAP, gamma=GAMMA, N=N_shaft),
        Sink(name="snk"),
    )


def _thermal_parts(comp, turb):
    comp_casing = ThermalMass(name="comp_casing", thermal_capacitance=200.0, T0=300.0)
    turb_casing = ThermalMass(name="turb_casing", thermal_capacitance=200.0, T0=300.0)
    shaft_mass = ThermalMass(name="shaft_mass", thermal_capacitance=50.0, T0=300.0)
    paths = (
        Convection(name="conv_comp_fluid", a=(comp, "out"), b=comp_casing, h=50.0, A=0.3),
        Convection(name="conv_turb_fluid", a=(turb, "out"), b=turb_casing, h=50.0, A=0.3),
        Conduction(name="cond_turb_shaft", a=turb_casing, b=shaft_mass, k=15.0, A=0.01, L=0.5),
        Conduction(name="cond_shaft_comp", a=shaft_mass, b=comp_casing, k=15.0, A=0.01, L=0.5),
        Convection(name="conv_turb_ambient", a=turb_casing, b=T_AMBIENT, h=10.0, A=1.0),
        Convection(name="conv_comp_ambient", a=comp_casing, b=T_AMBIENT, h=10.0, A=1.0),
    )
    return (comp_casing, turb_casing, shaft_mass), paths


def _wire_flow(network, src, comp, heater, turb, snk):
    for component in (src, comp, heater, turb, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", heater, "in")
    network.connect(heater, "out", turb, "in")
    network.connect(turb, "out", snk, "in")


def _build_via_add_heat_path():
    src, comp, heater, turb, snk = _flow_parts()
    masses, paths = _thermal_parts(comp, turb)
    network = Network(fluid=AIR)
    _wire_flow(network, src, comp, heater, turb, snk)
    for path in paths:
        network.add_heat_path(path)
    return network, comp, turb, masses, paths


def _build_by_hand():
    """The pre-add_heat_path() wiring, signs picked by hand — verbatim from
    tests/test_heat_transfer_integration.py."""
    src, comp, heater, turb, snk = _flow_parts()
    masses, paths = _thermal_parts(comp, turb)
    comp_casing, turb_casing, shaft_mass = masses
    (conv_comp_fluid, conv_turb_fluid, cond_turb_shaft,
     cond_shaft_comp, conv_turb_ambient, conv_comp_ambient) = paths

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
    _wire_flow(network, src, comp, heater, turb, snk)
    for component in (*masses, *paths):
        network.add_component(component)
    return network, comp, turb, masses, paths


def test_add_heat_path_reproduces_hand_wiring_exactly():
    # The whole point: auto-derived signs must match what a careful human
    # picks, and the resulting solve must be identical.
    auto_net, _, _, auto_masses, _ = _build_via_add_heat_path()
    hand_net, _, _, hand_masses, _ = _build_by_hand()

    for auto_mass, hand_mass in zip(auto_masses, hand_masses):
        auto_signs = {path.name: sign for path, sign in auto_mass.heat_sources}
        hand_signs = {path.name: sign for path, sign in hand_mass.heat_sources}
        assert auto_signs == hand_signs, f"{auto_mass.name} signs differ"

    auto = auto_net.solve(tol=1e-8, max_iter=300, damping=0.5, progress=False)
    hand = hand_net.solve(tol=1e-8, max_iter=300, damping=0.5, progress=False)
    for name in ("comp_casing.T", "turb_casing.T", "shaft_mass.T"):
        assert math.isclose(auto.params[name], hand.params[name], rel_tol=1e-9)


def test_add_heat_path_derives_mass_sign_from_which_endpoint_it_is():
    mass = ThermalMass(name="mass", thermal_capacitance=100.0, T0=300.0)
    network = Network(fluid=AIR)

    gains = Convection(name="gains", a=400.0, b=mass, h=1.0, A=1.0)  # mass is b
    loses = Convection(name="loses", a=mass, b=250.0, h=1.0, A=1.0)  # mass is a
    network.add_heat_path(gains)
    network.add_heat_path(loses)

    signs = {path.name: sign for path, sign in mass.heat_sources}
    assert signs == {"gains": 1.0, "loses": -1.0}


def test_add_heat_path_derives_flow_component_sign_from_which_endpoint_it_is():
    # Mirror image of the mass convention: heat_loss_watts() reports heat
    # *out of* the fluid, so a component is +1 as `a` (loses) and -1 as `b`.
    _, comp, _, turb, _ = _flow_parts()
    network = Network(fluid=AIR)

    losing = Convection(name="losing", a=(comp, "out"), b=300.0, h=1.0, A=1.0)
    gaining = Convection(name="gaining", a=900.0, b=(turb, "out"), h=1.0, A=1.0)
    network.add_heat_path(losing)
    network.add_heat_path(gaining)

    assert normalized_heat_paths(comp.heat_path) == [(losing, 1.0)]
    assert normalized_heat_paths(turb.heat_path) == [(gaining, -1.0)]


def test_add_heat_path_registers_thermal_mass_and_path_with_the_network():
    # A ThermalMass that never reaches add_component() contributes no
    # differential state, and the only symptom is a bare KeyError on its own
    # temperature from inside a residual call.
    mass = ThermalMass(name="mass", thermal_capacitance=100.0, T0=300.0)
    path = Convection(name="path", a=400.0, b=mass, h=1.0, A=1.0)
    network = Network(fluid=AIR)

    network.add_heat_path(path)

    assert mass in network.components
    assert path in network.components


def test_add_heat_path_is_idempotent():
    mass = ThermalMass(name="mass", thermal_capacitance=100.0, T0=300.0)
    path = Convection(name="path", a=400.0, b=mass, h=1.0, A=1.0)
    network = Network(fluid=AIR)

    network.add_heat_path(path)
    network.add_heat_path(path)

    assert mass.heat_sources == [(path, 1.0)]
    assert network.components.count(path) == 1


def test_add_heat_path_supports_several_paths_on_one_flow_component():
    # Previously impossible: heat_path was a single scalar attribute.
    _, comp, _, _, _ = _flow_parts()
    casing = ThermalMass(name="casing", thermal_capacitance=100.0, T0=300.0)
    network = Network(fluid=AIR)

    conv = Convection(name="conv", a=(comp, "out"), b=casing, h=50.0, A=0.3)
    rad = Radiation(name="rad", a=(comp, "out"), b=casing, emissivity=0.8, A=0.2)
    network.add_heat_path(conv)
    network.add_heat_path(rad)

    assert normalized_heat_paths(comp.heat_path) == [(conv, 1.0), (rad, 1.0)]
    assert {p.name for p, _ in casing.heat_sources} == {"conv", "rad"}


def test_heat_loss_watts_sums_every_attached_path_with_its_sign():
    class _FixedQ:
        def __init__(self, q):
            self.q = q

        def Q(self, state):
            return self.q

    assert heat_loss_watts(None, None) == 0.0
    assert heat_loss_watts(_FixedQ(120.0), None) == 120.0  # legacy scalar form
    assert heat_loss_watts([(_FixedQ(120.0), 1.0), (_FixedQ(50.0), -1.0)], None) == 70.0
    assert heat_loss_watts([_FixedQ(120.0), _FixedQ(50.0)], None) == 170.0  # bare, +1 each


def test_add_heat_path_rejects_a_bad_port_name_at_wiring_time():
    # Previously a bare KeyError raised deep inside _temperature_of during a
    # residual evaluation, long after the mistake was made.
    _, comp, _, _, _ = _flow_parts()
    casing = ThermalMass(name="casing", thermal_capacitance=100.0, T0=300.0)
    network = Network(fluid=AIR)

    path = Convection(name="path", a=(comp, "outlet"), b=casing, h=50.0, A=0.3)
    with pytest.raises(NetworkTopologyError, match="no port 'outlet'"):
        network.add_heat_path(path)


def test_add_heat_path_rejects_a_component_with_no_heat_path_attribute():
    _, _, heater, _, _ = _flow_parts()  # Pipe takes a scalar heat_loss instead
    casing = ThermalMass(name="casing", thermal_capacitance=100.0, T0=300.0)
    network = Network(fluid=AIR)

    path = Convection(name="path", a=(heater, "out"), b=casing, h=50.0, A=0.3)
    with pytest.raises(NetworkTopologyError, match="no heat_path attribute"):
        network.add_heat_path(path)


def test_add_heat_path_rejects_something_that_is_not_a_heat_path():
    network = Network(fluid=AIR)
    with pytest.raises(NetworkTopologyError, match="endpoints"):
        network.add_heat_path(ThermalMass(name="not_a_path", thermal_capacitance=1.0, T0=300.0))


def test_thermal_masses_wired_by_add_heat_path_reach_zero_net_heat_at_steady_state():
    # The physical guard: at steady state every mass's signed heat sum must
    # be zero. A flipped sign still solves, but balances the wrong equation.
    network, _, _, masses, _ = _build_via_add_heat_path()

    result = network.solve(tol=1e-8, max_iter=300, damping=0.5, progress=False)
    state = result.state()

    for mass in masses:
        assert math.isclose(mass._net_heat(state), 0.0, abs_tol=1e-6)
