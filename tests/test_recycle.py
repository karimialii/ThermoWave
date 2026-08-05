import pytest

pytest.importorskip("CoolProp")

from thermowave.components.pump import Pump  # noqa: E402
from thermowave.components.recycle import Recycle  # noqa: E402
from thermowave.components.simple_condenser import SimpleCondenser  # noqa: E402
from thermowave.components.simple_evaporator import SimpleEvaporator  # noqa: E402
from thermowave.components.steam_turbine import SteamTurbine  # noqa: E402
from thermowave.core.network import Network  # noqa: E402
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


def test_rejects_nonpositive_mdot():
    with pytest.raises(ValueError, match="mdot"):
        Recycle(name="rc", mdot=0.0)
    with pytest.raises(ValueError, match="mdot"):
        Recycle(name="rc", mdot=-1.0)


def test_ports_and_fixed_mdot():
    rc = Recycle(name="rc", mdot=2.5)
    assert rc.ports() == {"in": "rc.in", "out": "rc.out"}
    assert rc.fixed_node_mdot() == {"rc.out": 2.5}


def test_residuals_pass_through_hand_calc():
    rc = Recycle(name="rc", mdot=1.0)
    P_in, h_in = 1.0e5, 5.0e5
    P_out, h_out = 1.2e5, 5.5e5
    state = _FakeState(
        fluid=WATER,
        node_values={"rc.in": (P_in, h_in), "rc.out": (P_out, h_out)},
        mdots={},
    )
    residuals = rc.residuals(state)
    assert residuals == pytest.approx([P_out - P_in, h_out - h_in])


def test_residuals_zero_when_out_equals_in():
    rc = Recycle(name="rc", mdot=1.0)
    P, h = 1.0e5, 5.0e5
    state = _FakeState(
        fluid=WATER, node_values={"rc.in": (P, h), "rc.out": (P, h)}, mdots={}
    )
    assert rc.residuals(state) == pytest.approx([0.0, 0.0])


def test_closed_rankine_loop_with_no_source_or_sink_converges():
    net = Network(fluid=WATER)
    mdot = 1.0
    pump = Pump(name="pump", P_out=5.0e6, eta=0.75)
    evap = SimpleEvaporator(name="evap", superheat=30.0)
    turb = SteamTurbine(name="turb", P_out=1.0e4, eta_s=0.85)
    cond = SimpleCondenser(name="cond", subcool=2.0)
    recycle = Recycle(name="recycle", mdot=mdot)

    for c in (pump, evap, turb, cond, recycle):
        net.add_component(c)
    net.connect(pump, "out", evap, "in")
    net.connect(evap, "out", turb, "in")
    net.connect(turb, "out", cond, "in")
    net.connect(cond, "out", recycle, "in")
    net.connect(recycle, "out", pump, "in")

    result = net.solve(progress=False)

    assert result.converged
    for node in ("pump.out", "evap.out", "turb.out", "cond.out", "recycle.out"):
        assert result.node_mdot[node] == pytest.approx(mdot)

    state = result.state()
    turb_power = turb.report_metrics(state)["power [W]"]
    pump_power = pump.report_metrics(state)["power [W]"]
    q_evap = evap.report_metrics(state)["power [W]"]
    q_cond = -cond.report_metrics(state)["power [W]"]  # negative: heat rejected

    # Closed-cycle energy balance: heat in - heat out == net shaft work out.
    assert (q_evap - q_cond) == pytest.approx(turb_power - pump_power, rel=1e-6)


class _FakeCompositionFluid:
    """Minimal stand-in exposing just the mass_fractions() duck-type
    Recycle.recycle_fluid_delta() relies on -- no need for a real
    CoolProp/Cantera fluid to unit-test the tear-stream bookkeeping itself."""

    def __init__(self, Y):
        self._Y = Y

    def mass_fractions(self):
        return dict(self._Y)


class _NoCompositionFluid:
    """A fluid exposing no mass_fractions() at all -- the "nothing
    comparable" case recycle_fluid_delta() must fall back on."""


class _FakeFluidState:
    def __init__(self, fluid_at_map):
        self._map = fluid_at_map

    def fluid_at(self, name):
        return self._map[name]


def test_fluid_seed_empty_without_guess():
    rc = Recycle(name="rc", mdot=1.0)
    assert rc.fluid_seed() == {}


def test_fluid_seed_returns_guess_at_outlet_node():
    guess = _FakeCompositionFluid({"N2": 1.0})
    rc = Recycle(name="rc", mdot=1.0, fluid_guess=guess)
    assert rc.fluid_seed() == {"rc.out": guess}


def test_recycle_fluid_delta_none_without_guess():
    rc = Recycle(name="rc", mdot=1.0)
    state = _FakeFluidState({"rc.in": _FakeCompositionFluid({"N2": 1.0})})
    assert rc.recycle_fluid_delta(state) is None


def test_recycle_fluid_delta_none_when_not_comparable():
    rc = Recycle(name="rc", mdot=1.0, fluid_guess=_NoCompositionFluid())
    state = _FakeFluidState({"rc.in": _NoCompositionFluid()})
    assert rc.recycle_fluid_delta(state) is None


def test_recycle_fluid_delta_is_max_mass_fraction_gap():
    guess = _FakeCompositionFluid({"N2": 0.79, "O2": 0.21})
    actual = _FakeCompositionFluid({"N2": 0.70, "O2": 0.21, "CO2": 0.09})
    rc = Recycle(name="rc", mdot=1.0, fluid_guess=guess)
    state = _FakeFluidState({"rc.in": actual})
    assert rc.recycle_fluid_delta(state) == pytest.approx(0.09)


def test_recycle_fluid_delta_zero_when_guess_already_matches():
    fluid = _FakeCompositionFluid({"N2": 0.79, "O2": 0.21})
    rc = Recycle(name="rc", mdot=1.0, fluid_guess=fluid)
    state = _FakeFluidState({"rc.in": fluid})
    assert rc.recycle_fluid_delta(state) == pytest.approx(0.0)


def test_update_fluid_guess_replaces_guess_with_inlet_fluid():
    guess = _FakeCompositionFluid({"N2": 1.0})
    actual = _FakeCompositionFluid({"CO2": 1.0})
    rc = Recycle(name="rc", mdot=1.0, fluid_guess=guess)
    state = _FakeFluidState({"rc.in": actual})
    rc.update_fluid_guess(state)
    assert rc._fluid_guess is actual


def test_update_fluid_guess_is_noop_without_guess():
    rc = Recycle(name="rc", mdot=1.0)
    state = _FakeFluidState({"rc.in": _FakeCompositionFluid({"N2": 1.0})})
    rc.update_fluid_guess(state)
    assert rc._fluid_guess is None
